// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <algorithm>
#include "BQPDSolver.hpp"
#include "optimization/Direction.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "fortran_interface.h"

#define WSC FC_GLOBAL(wsc, WSC)
#define ALPHAC FC_GLOBAL(alphac, ALPHAC)
#define BQPD FC_GLOBAL(bqpd, BQPD)
#define hessian_vector_product FC_GLOBAL(gdotx, GDOTX)

extern "C" {
   void hessian_vector_product([[maybe_unused]] int *n, const double x[], const double ws[], const int lws[], double v[]);

   // fortran common block used in bqpd/bqpd.f
   extern struct {
      int kk, ll, kkk, lll, mxws, mxlws;
   } WSC;

   // fortran common for inertia correction in wdotd
   extern struct {
      double alpha;
   } ALPHAC;

   extern void
   BQPD(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
         double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
         int* info, int* iprint, int* nout);
}

namespace uno {
   #define BIG 1e30

   // preallocate a bunch of stuff
   BQPDSolver::BQPDSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
         size_t number_hessian_nonzeros, BQPDProblemType problem_type, const Options& options):
         QPSolver(), number_hessian_nonzeros(number_hessian_nonzeros),
         lb(number_variables + number_constraints),
         ub(number_variables + number_constraints),
         jacobian(number_jacobian_nonzeros + number_objective_gradient_nonzeros), // Jacobian + objective gradient
         jacobian_sparsity(number_jacobian_nonzeros + number_objective_gradient_nonzeros + number_constraints + 3),
         kmax(problem_type == BQPDProblemType::QP ? options.get_int("BQPD_kmax") : 0), alp(static_cast<size_t>(this->mlp)),
         lp(static_cast<size_t>(this->mlp)),
         active_set(number_variables + number_constraints),
         w(number_variables + number_constraints), gradient_solution(number_variables), residuals(number_variables + number_constraints),
         e(number_variables + number_constraints),
         size_hessian_sparsity(problem_type == BQPDProblemType::QP ? number_hessian_nonzeros + number_variables + 3 : 0),
         size_hessian_workspace(number_hessian_nonzeros + static_cast<size_t>(this->kmax * (this->kmax + 9) / 2) + 2 * number_variables +
                                number_constraints + this->mxwk0),
         size_hessian_sparsity_workspace(this->size_hessian_sparsity + static_cast<size_t>(this->kmax) + this->mxiwk0),
         hessian_values(this->size_hessian_workspace),
         hessian_sparsity(this->size_hessian_sparsity_workspace),
         current_hessian_indices(number_variables),
         print_subproblem(options.get_bool("print_subproblem")) {
      // default active set
      for (size_t variable_index: Range(number_variables + number_constraints)) {
         this->active_set[variable_index] = static_cast<int>(variable_index) + this->fortran_shift;
      }
   }

   void BQPDSolver::solve_QP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<size_t, double>& hessian, const Vector<double>& initial_point,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->save_hessian_to_local_format(hessian);
      }
      if (this->print_subproblem) {
         DEBUG << "QP:\n";
         DEBUG << "Hessian: " << hessian;
      }
      this->solve_subproblem(number_variables, number_constraints, variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
            constraints_upper_bounds, linear_objective, constraint_jacobian, initial_point, direction, warmstart_information);
   }

   void BQPDSolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      if (this->print_subproblem) {
         DEBUG << "LP:\n";
      }
      this->solve_subproblem(number_variables, number_constraints, variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds,
            constraints_upper_bounds, linear_objective, constraint_jacobian, initial_point, direction, warmstart_information);
   }

   void BQPDSolver::solve_subproblem(size_t number_variables, size_t number_constraints, const std::vector<double>& variables_lower_bounds,
         const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
         const std::vector<double>& constraints_upper_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) {
      // initialize wsc_ common block (Hessian & workspace for BQPD)
      // setting the common block here ensures that several instances of BQPD can run simultaneously
      WSC.kk = static_cast<int>(this->number_hessian_nonzeros);
      WSC.ll = static_cast<int>(this->size_hessian_sparsity);
      WSC.mxws = static_cast<int>(this->size_hessian_workspace);
      WSC.mxlws = static_cast<int>(this->size_hessian_sparsity_workspace);
      ALPHAC.alpha = 0.; // inertia control

      if (this->print_subproblem) {
         DEBUG << "objective gradient: " << linear_objective;
         for (size_t constraint_index: Range(number_constraints)) {
            DEBUG << "gradient c" << constraint_index << ": " << constraint_jacobian[constraint_index];
         }
         for (size_t variable_index: Range(number_variables)) {
            DEBUG << "d" << variable_index << " in [" << variables_lower_bounds[variable_index] << ", " << variables_upper_bounds[variable_index] << "]\n";
         }
         for (size_t constraint_index: Range(number_constraints)) {
            DEBUG << "linearized c" << constraint_index << " in [" << constraints_lower_bounds[constraint_index] << ", " << constraints_upper_bounds[constraint_index] << "]\n";
         }
      }

      // Jacobian (objective and constraints)
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->save_gradients_to_local_format(number_constraints, linear_objective, constraint_jacobian);
      }

      // set variable bounds
      if (warmstart_information.variable_bounds_changed) {
         for (size_t variable_index: Range(number_variables)) {
            this->lb[variable_index] = (variables_lower_bounds[variable_index] == -INF<double>) ? -BIG : variables_lower_bounds[variable_index];
            this->ub[variable_index] = (variables_upper_bounds[variable_index] == INF<double>) ? BIG : variables_upper_bounds[variable_index];
         }
      }
      // set constraint bounds
      if (warmstart_information.constraint_bounds_changed) {
         for (size_t constraint_index: Range(number_constraints)) {
            this->lb[number_variables + constraint_index] = (constraints_lower_bounds[constraint_index] == -INF<double>) ? -BIG : constraints_lower_bounds[constraint_index];
            this->ub[number_variables + constraint_index] = (constraints_upper_bounds[constraint_index] == INF<double>) ? BIG : constraints_upper_bounds[constraint_index];
         }
      }

      direction.primals = initial_point;
      const int n = static_cast<int>(number_variables);
      const int m = static_cast<int>(number_constraints);

      const BQPDMode mode = this->determine_mode(warmstart_information);
      const int mode_integer = static_cast<int>(mode);

      // solve the LP/QP
      BQPD(&n, &m, &this->k, &this->kmax, this->jacobian.data(), this->jacobian_sparsity.data(), direction.primals.data(), this->lb.data(),
            this->ub.data(), &direction.subproblem_objective, &this->fmin, this->gradient_solution.data(), this->residuals.data(), this->w.data(),
            this->e.data(), this->active_set.data(), this->alp.data(), this->lp.data(), &this->mlp, &this->peq_solution, this->hessian_values.data(),
            this->hessian_sparsity.data(), &mode_integer, &this->ifail, this->info.data(), &this->iprint, &this->nout);
      const BQPDStatus bqpd_status = BQPDSolver::bqpd_status_from_int(this->ifail);
      direction.status = BQPDSolver::status_from_bqpd_status(bqpd_status);
      this->number_calls++;

      // project solution into bounds
      for (size_t variable_index: Range(number_variables)) {
         direction.primals[variable_index] = std::min(std::max(direction.primals[variable_index], variables_lower_bounds[variable_index]),
               variables_upper_bounds[variable_index]);
      }
      this->set_multipliers(number_variables, direction.multipliers);
   }

   BQPDMode BQPDSolver::determine_mode(const WarmstartInformation& warmstart_information) const {
      BQPDMode mode = (this->number_calls == 0) ? BQPDMode::ACTIVE_SET_EQUALITIES : BQPDMode::USER_DEFINED;
      // if problem structure changed, use cold start
      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         mode = BQPDMode::ACTIVE_SET_EQUALITIES;
      }
      // if only the variable bounds changed, reuse the active set estimate and the Jacobian information
      else if (warmstart_information.variable_bounds_changed && not warmstart_information.objective_changed &&
               not warmstart_information.constraints_changed && not warmstart_information.constraint_bounds_changed) {
         mode = BQPDMode::UNCHANGED_ACTIVE_SET_AND_JACOBIAN;
      }
      return mode;
   }

   // save Hessian (in arbitrary format) to a "weak" CSC format: compressed columns but row indices are not sorted, nor unique
   void BQPDSolver::save_hessian_to_local_format(const SymmetricMatrix<size_t, double>& hessian) {
      const size_t header_size = 1;
      // pointers withing the single array
      int* row_indices = &this->hessian_sparsity[header_size];
      int* column_starts = &this->hessian_sparsity[header_size + hessian.number_nonzeros()];
      // header
      this->hessian_sparsity[0] = static_cast<int>(hessian.number_nonzeros() + 1);
      // count the elements in each column
      for (size_t column_index: Range(hessian.dimension() + 1)) {
         column_starts[column_index] = 0;
      }
      for (const auto [row_index, column_index, element]: hessian) {
         column_starts[column_index + 1]++;
      }
      // carry over the column starts
      for (size_t column_index: Range(1, hessian.dimension() + 1)) {
         column_starts[column_index] += column_starts[column_index - 1];
         column_starts[column_index - 1] += this->fortran_shift;
      }
      column_starts[hessian.dimension()] += this->fortran_shift;
      // copy the entries
      //std::vector<int> current_indices(hessian.dimension());
      this->current_hessian_indices.fill(0);
      for (const auto [row_index, column_index, element]: hessian) {
         const size_t index = static_cast<size_t>(column_starts[column_index] + this->current_hessian_indices[column_index] - this->fortran_shift);
         assert(index <= static_cast<size_t>(column_starts[column_index + 1]) &&
                "BQPD: error in converting the Hessian matrix to the local format. Try setting the sparse format to CSC");
         this->hessian_values[index] = element;
         row_indices[index] = static_cast<int>(row_index) + this->fortran_shift;
         this->current_hessian_indices[column_index]++;
      }
   }

   void BQPDSolver::save_gradients_to_local_format(size_t number_constraints, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian) {
      size_t current_index = 0;
      for (const auto [variable_index, derivative]: linear_objective) {
         this->jacobian[current_index] = derivative;
         this->jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
         current_index++;
      }
      for (size_t constraint_index: Range(number_constraints)) {
         for (const auto [variable_index, derivative]: constraint_jacobian[constraint_index]) {
            this->jacobian[current_index] = derivative;
            this->jacobian_sparsity[current_index + 1] = static_cast<int>(variable_index) + this->fortran_shift;
            current_index++;
         }
      }
      current_index++;
      this->jacobian_sparsity[0] = static_cast<int>(current_index);
      // header
      size_t size = 1;
      this->jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
      size += linear_objective.size();
      this->jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
      for (size_t constraint_index: Range(number_constraints)) {
         size += constraint_jacobian[constraint_index].size();
         this->jacobian_sparsity[current_index] = static_cast<int>(size);
         current_index++;
      }
   }

   void BQPDSolver::set_multipliers(size_t number_variables, Multipliers& direction_multipliers) {
      direction_multipliers.reset();

      // active constraints
      for (size_t active_constraint_index: Range(number_variables - static_cast<size_t>(this->k))) {
         const size_t index = static_cast<size_t>(std::abs(this->active_set[active_constraint_index]) - this->fortran_shift);

         if (index < number_variables) {
            // bound constraint
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction_multipliers.lower_bounds[index] = this->residuals[index];
            }
            else { // upper bound active */
               direction_multipliers.upper_bounds[index] = -this->residuals[index];
            }
         }
         else {
            // general constraint
            size_t constraint_index = index - number_variables;
            if (0 <= this->active_set[active_constraint_index]) { // lower bound active
               direction_multipliers.constraints[constraint_index] = this->residuals[index];
            }
            else { // upper bound active
               direction_multipliers.constraints[constraint_index] = -this->residuals[index];
            }
         }
      }
   }

   BQPDStatus BQPDSolver::bqpd_status_from_int(int ifail) {
      assert(0 <= ifail && ifail <= 9 && "BQPDSolver.bqpd_status_from_int: ifail does not belong to [0, 9]");
      return static_cast<BQPDStatus>(ifail);
   }

   SubproblemStatus BQPDSolver::status_from_bqpd_status(BQPDStatus bqpd_status) {
      switch (bqpd_status) {
         case BQPDStatus::OPTIMAL:
            return SubproblemStatus::OPTIMAL;
         case BQPDStatus::UNBOUNDED_PROBLEM:
            return SubproblemStatus::UNBOUNDED_PROBLEM;
         case BQPDStatus::BOUND_INCONSISTENCY:
            DEBUG << "BQPD error: bound inconsistency\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::INFEASIBLE:
            return SubproblemStatus::INFEASIBLE;
            // errors
         case BQPDStatus::INCORRECT_PARAMETER:
            DEBUG << "BQPD error: incorrect parameter\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::LP_INSUFFICIENT_SPACE:
            DEBUG << "BQPD error: LP insufficient space\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::HESSIAN_INSUFFICIENT_SPACE:
            DEBUG << "BQPD kmax too small, continue anyway\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::SPARSE_INSUFFICIENT_SPACE:
            DEBUG << "BQPD error: sparse insufficient space\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::MAX_RESTARTS_REACHED:
            DEBUG << "BQPD max restarts reached\n";
            return SubproblemStatus::ERROR;
         case BQPDStatus::UNDEFINED:
            DEBUG << "BQPD error: undefined\n";
            return SubproblemStatus::ERROR;
      }
      throw std::invalid_argument("The BQPD ifail is not consistent with the Uno status values");
   }
} // namespace

void hessian_vector_product(int *n, const double x[], const double ws[], const int lws[], double v[]) {
   assert(n != nullptr && "BQPDSolver::hessian_vector_product: the dimension n passed by pointer is NULL");

   for (size_t i = 0; i < static_cast<size_t>(*n); i++) {
      v[i] = 0.;
   }

   int footer_start = lws[0];
   for (int i = 0; i < *n; i++) {
      for (int k = lws[footer_start + i]; k < lws[footer_start + i + 1]; k++) {
         int j = lws[k] - 1;
         v[i] += ws[k - 1] * x[j];
         if (j != i) {
            // off-diagonal term
            v[j] += ws[k - 1] * x[i];
         }
      }
   }
}