#include <cassert>
#include "HiGHSSolver.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   HiGHSSolver::HiGHSSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, size_t /*number_hessian_nonzeros*/, const Options& options):
         LPSolver(),
         constraints(number_constraints),
         linear_objective(number_objective_gradient_nonzeros),
         constraint_jacobian(number_constraints, number_variables),
         print_subproblem(options.get_bool("print_subproblem")) {
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(number_variables);
      this->model.lp_.col_upper_.resize(number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(number_constraints);
      this->model.lp_.row_upper_.resize(number_constraints);
      // constraint matrix in CSR format (each row is a constraint gradient)
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
      this->model.lp_.a_matrix_.value_.reserve(number_jacobian_nonzeros);
      this->model.lp_.a_matrix_.index_.reserve(number_jacobian_nonzeros);
      this->model.lp_.a_matrix_.start_.reserve(number_variables + 1);

      this->highs_solver.setOptionValue("output_flag", "false");
   }

   void HiGHSSolver::solve_LP(Statistics& /*statistics*/, LagrangeNewtonSubproblem& subproblem, const Vector<double>& /*initial_point*/,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      if (this->print_subproblem) {
         DEBUG << "LP:\n";
      }
      this->set_up_subproblem(subproblem, warmstart_information);
      this->solve_subproblem(subproblem, direction);
   }

   void HiGHSSolver::set_up_subproblem(LagrangeNewtonSubproblem& subproblem, const WarmstartInformation& warmstart_information) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(subproblem.number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(subproblem.number_constraints);

      // function evaluations
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->linear_objective);
         // reset the linear part of the objective
         for (size_t variable_index: Range(subproblem.number_variables)) {
            this->model.lp_.col_cost_[variable_index] = 0.;
         }
         // copy the sparse vector into the dense vector
         for (const auto [variable_index, derivative]: this->linear_objective) {
            this->model.lp_.col_cost_[variable_index] = derivative;
         }
      }
      if (warmstart_information.constraints_changed) {
         // constraints
         subproblem.evaluate_constraints(this->constraints);

         // constraint Jacobian
         subproblem.evaluate_constraint_jacobian(this->constraint_jacobian);
         this->model.lp_.a_matrix_.value_.clear();
         this->model.lp_.a_matrix_.index_.clear();
         this->model.lp_.a_matrix_.start_.clear();

         size_t number_nonzeros = 0;
         this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
         for (size_t constraint_index: Range(subproblem.number_constraints)) {
            for (const auto [variable_index, value]: this->constraint_jacobian[constraint_index]) {
               this->model.lp_.a_matrix_.value_.emplace_back(value);
               this->model.lp_.a_matrix_.index_.emplace_back(variable_index);
               number_nonzeros++;
            }
            this->model.lp_.a_matrix_.start_.emplace_back(number_nonzeros);
         }
      }

      // variable bounds
      if (warmstart_information.variable_bounds_changed) {
         subproblem.set_variable_bounds(this->model.lp_.col_lower_, this->model.lp_.col_upper_);
      }

      // constraint bounds
      if (warmstart_information.constraint_bounds_changed || warmstart_information.constraints_changed) {
         subproblem.set_constraint_bounds(this->constraints, this->model.lp_.row_lower_, this->model.lp_.row_upper_);
      }

      if (this->print_subproblem) {
         DEBUG << "Linear objective part: "; print_vector(DEBUG, view(this->model.lp_.col_cost_, 0, subproblem.number_variables));
         DEBUG << "Jacobian:\n";
         DEBUG << "J = "; print_vector(DEBUG, this->model.lp_.a_matrix_.value_);
         DEBUG << "with column start: "; print_vector(DEBUG, this->model.lp_.a_matrix_.start_);
         DEBUG << "and row index: "; print_vector(DEBUG, this->model.lp_.a_matrix_.index_);
         for (size_t variable_index = 0; variable_index < subproblem.number_variables; variable_index++) {
            DEBUG << "d" << variable_index << " in [" << this->model.lp_.col_lower_[variable_index] << ", " <<
               this->model.lp_.col_upper_[variable_index] << "]\n";
         }
         for (size_t constraint_index = 0; constraint_index < subproblem.number_constraints; constraint_index++) {
            DEBUG << "linearized c" << constraint_index << " in [" << this->model.lp_.row_lower_[constraint_index] << ", " <<
               this->model.lp_.row_upper_[constraint_index]<< "]\n";
         }
      }
   }

   void HiGHSSolver::solve_subproblem(LagrangeNewtonSubproblem& subproblem, Direction& direction) {
      // solve the LP
      HighsStatus return_status = this->highs_solver.passModel(this->model);
      assert(return_status == HighsStatus::kOk);

      return_status = this->highs_solver.run(); // solve
      DEBUG << "HiGHS status: " << static_cast<int>(return_status) << '\n';

      // if HiGHS could not optimize (e.g. because of indefinite Hessian), return an error
      if (return_status == HighsStatus::kError) {
         direction.status = SubproblemStatus::ERROR;
         return;
      }
      HighsModelStatus model_status = highs_solver.getModelStatus();
      DEBUG << "HiGHS model status: " << static_cast<int>(model_status) << '\n';

      if (model_status == HighsModelStatus::kInfeasible) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      else if (model_status == HighsModelStatus::kUnbounded) {
         direction.status = SubproblemStatus::UNBOUNDED_PROBLEM;
         return;
      }
      
      direction.status = SubproblemStatus::OPTIMAL;
      const HighsSolution& solution = this->highs_solver.getSolution();
      // read the primal solution and bound dual solution
      for (size_t variable_index = 0; variable_index < subproblem.number_variables; variable_index++) {
         direction.primals[variable_index] = solution.col_value[variable_index];
         const double bound_multiplier = solution.col_dual[variable_index];
         if (0. < bound_multiplier) {
            direction.multipliers.lower_bounds[variable_index] = bound_multiplier;
         }
         else {
            direction.multipliers.upper_bounds[variable_index] = bound_multiplier;
         }
      }
      // read the dual solution
      for (size_t constraint_index = 0; constraint_index < subproblem.number_constraints; constraint_index++) {
         direction.multipliers.constraints[constraint_index] = solution.row_dual[constraint_index];
      }
      const HighsInfo& info = this->highs_solver.getInfo();
      direction.subproblem_objective = info.objective_function_value;
   }
} // namespace