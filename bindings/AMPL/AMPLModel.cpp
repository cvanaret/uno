// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "AMPLModel.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"
#include "tools/Infinity.hpp"
#include "symbolic/Concatenation.hpp"

namespace uno {
   ASL* generate_asl(std::string file_name) {
      ASL* asl = ASL_alloc(ASL_read_pfgh);
      FILE* nl = jac0dim_ASL(asl, file_name.data(), static_cast<int>(file_name.size()));
      // indices start at 0
      asl->i.Fortran_ = 0;

      int n_discrete = asl->i.nbv_ + asl->i.niv_ + asl->i.nlvbi_ + asl->i.nlvci_ + asl->i.nlvoi_;
      if (0 < n_discrete) {
         WARNING << "Ignoring integrality of " << n_discrete << " variables.\n";
         asl->i.need_nl_ = 0;
      }

      // preallocate initial primal and dual solutions
      asl->i.X0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(asl->i.n_var_)));
      asl->i.pi0_ = static_cast<double*>(M1zapalloc_ASL(&asl->i, sizeof(double) * static_cast<size_t>(asl->i.n_con_)));

      // read the file_name.nl file
      pfgh_read_ASL(asl, nl, ASL_findgroups);
      return asl;
   }

   // generate the ASL object and call the private constructor
   AMPLModel::AMPLModel(const std::string& file_name) : AMPLModel(file_name, generate_asl(file_name)) {
   }

   AMPLModel::AMPLModel(const std::string& file_name, ASL* asl) :
         Model(file_name, static_cast<size_t>(asl->i.n_var_), static_cast<size_t>(asl->i.n_con_), (asl->i.objtype_[0] == 1) ? -1. : 1.),
         asl(asl),
         // allocate vectors
         asl_gradient(this->number_variables),
         variable_lower_bounds(this->number_variables),
         variable_upper_bounds(this->number_variables),
         constraint_lower_bounds(this->number_constraints),
         constraint_upper_bounds(this->number_constraints),
         variable_status(this->number_variables),
         constraint_type(this->number_constraints),
         constraint_status(this->number_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints),
         lower_bounded_variables_collection(this->lower_bounded_variables),
         upper_bounded_variables_collection(this->upper_bounded_variables),
         single_lower_bounded_variables_collection(this->single_lower_bounded_variables),
         single_upper_bounded_variables_collection(this->single_upper_bounded_variables) {
      // evaluate the constraint Jacobian in sparse mode
      this->asl->i.congrd_mode = 1;

      // variables
      this->lower_bounded_variables.reserve(this->number_variables);
      this->upper_bounded_variables.reserve(this->number_variables);
      this->single_lower_bounded_variables.reserve(this->number_variables);
      this->single_upper_bounded_variables.reserve(this->number_variables);
      this->generate_variables();

      // constraints
      this->equality_constraints.reserve(this->number_constraints);
      this->inequality_constraints.reserve(this->number_constraints);
      this->linear_constraints.reserve(this->number_constraints);
      this->generate_constraints();

      // compute number of nonzeros in the Lagrangian Hessian
      this->set_number_hessian_nonzeros();
   }

   AMPLModel::~AMPLModel() {
      ASL_free(&this->asl);
   }

   void AMPLModel::generate_variables() {
      for (size_t variable_index: Range(this->number_variables)) {
         this->variable_lower_bounds[variable_index] = (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2*variable_index] : -INF<double>;
         this->variable_upper_bounds[variable_index] = (this->asl->i.LUv_ != nullptr) ? this->asl->i.LUv_[2*variable_index + 1] : INF<double>;
         if (this->variable_lower_bounds[variable_index] == this->variable_upper_bounds[variable_index]) {
            WARNING << "Variable x" << variable_index << " has identical bounds\n";
         }
      }
      AMPLModel::determine_bounds_types(this->variable_lower_bounds, this->variable_upper_bounds, this->variable_status);
      // figure out the bounded variables
      for (size_t variable_index: Range(this->number_variables)) {
         const BoundType status = this->get_variable_bound_type(variable_index);
         if (status == BOUNDED_LOWER || status == BOUNDED_BOTH_SIDES) {
            this->lower_bounded_variables.push_back(variable_index);
            if (status == BOUNDED_LOWER) {
               this->single_lower_bounded_variables.push_back(variable_index);
            }
         }
         if (status == BOUNDED_UPPER || status == BOUNDED_BOTH_SIDES) {
            this->upper_bounded_variables.push_back(variable_index);
            if (status == BOUNDED_UPPER) {
               this->single_upper_bounded_variables.push_back(variable_index);
            }
         }
      }
   }

   double AMPLModel::evaluate_objective(const Vector<double>& x) const {
      int error_flag = 0;
      double result = this->objective_sign * (*(this->asl)->p.Objval)(this->asl, 0, const_cast<double*>(x.data()), &error_flag);
      if (0 < error_flag) {
         throw FunctionEvaluationError();
      }
      return result;
   }

   // sparse gradient
   void AMPLModel::evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const {
      int error_flag = 0;
      // prevent ASL to crash by catching all evaluation errors
      Jmp_buf err_jmp_uno;
      asl->i.err_jmp_ = &err_jmp_uno;
      asl->i.err_jmp1_ = &err_jmp_uno;
      if (setjmp(err_jmp_uno.jb)) {
         error_flag = 1;
      }
      // evaluate the ASL gradient (always in a dense vector)
      (*(this->asl)->p.Objgrd)(this->asl, 0, const_cast<double*>(x.data()), const_cast<double*>(this->asl_gradient.data()), &error_flag);
      if (0 < error_flag) {
         throw GradientEvaluationError();
      }

      // create the Uno sparse vector
      ograd* asl_variables_tmp = this->asl->i.Ograd_[0];
      while (asl_variables_tmp != nullptr) {
         const size_t variable_index = static_cast<size_t>(asl_variables_tmp->varno);
         // scale by the objective sign
         const double partial_derivative = this->objective_sign*this->asl_gradient[variable_index];
         gradient.insert(variable_index, partial_derivative);
         asl_variables_tmp = asl_variables_tmp->next;
      }
   }

   /*
   double AMPLModel::evaluate_constraint(int j, const std::vector<double>& x) const {
      int error_flag = 0;
      double result = (*(this->asl)->p.Conival)(this->asl_, j, const_cast<double*>(x.data()), &error_flag);
      if (0 < error_flag) {
         throw FunctionNumericalError();
      }
      return result;
   }
   */

   void AMPLModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      int error_flag = 0;
      (*(this->asl)->p.Conval)(this->asl, const_cast<double*>(x.data()), constraints.data(), &error_flag);
      if (0 < error_flag) {
         throw FunctionEvaluationError();
      }
   }

   // sparse gradient
   void AMPLModel::evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
      // compute the AMPL sparse gradient
      int error_flag = 0;
      (*(this->asl)->p.Congrd)(this->asl, static_cast<int>(constraint_index), const_cast<double*>(x.data()), const_cast<double*>(this->asl_gradient.data()),
            &error_flag);
      if (0 < error_flag) {
         throw GradientEvaluationError();
      }

      // construct the Uno sparse vector
      gradient.clear();
      cgrad* asl_variables_tmp = this->asl->i.Cgrad_[constraint_index];
      size_t sparse_asl_index = 0;
      while (asl_variables_tmp != nullptr) {
         const size_t variable_index = static_cast<size_t>(asl_variables_tmp->varno);
         gradient.insert(variable_index, this->asl_gradient[sparse_asl_index]);
         asl_variables_tmp = asl_variables_tmp->next;
         sparse_asl_index++;
      }
   }

   void AMPLModel::evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
      for (size_t constraint_index: Range(this->number_constraints)) {
         constraint_jacobian[constraint_index].clear();
         this->evaluate_constraint_gradient(x, constraint_index, constraint_jacobian[constraint_index]);
      }
   }

   void AMPLModel::set_number_hessian_nonzeros() {
      // compute the maximum number of nonzero elements, provided that all multipliers are non-zero
      // int (*Sphset) (ASL*, SputInfo**, int nobj, int ow, int y, int uptri);
      const int objective_number = -1;
      const int upper_triangular = 1;
      this->number_asl_hessian_nonzeros = static_cast<size_t>((*(this->asl)->p.Sphset)(this->asl, nullptr, objective_number, 1, 1, upper_triangular));
      this->asl_hessian.reserve(this->number_asl_hessian_nonzeros);

      // use Lagrangian scale: in AMPL, the Lagrangian is f + lambda.g, while Uno uses f - lambda.g
      int error_flag{};
      lagscale_ASL(this->asl, -1., &error_flag);
   }

   size_t AMPLModel::number_objective_gradient_nonzeros() const {
      return static_cast<size_t>(this->asl->i.nzo_);
   }

   size_t AMPLModel::number_jacobian_nonzeros() const {
      return static_cast<size_t>(this->asl->i.nzc_);
   }

   size_t AMPLModel::number_hessian_nonzeros() const {
      return this->number_asl_hessian_nonzeros;
   }

   const Collection<size_t>& AMPLModel::get_equality_constraints() const {
      return this->equality_constraints_collection;
   }

   const Collection<size_t>& AMPLModel::get_inequality_constraints() const {
      return this->inequality_constraints_collection;
   }

   const std::vector<size_t>& AMPLModel::get_linear_constraints() const {
      return this->linear_constraints;
   }

   const SparseVector<size_t>& AMPLModel::get_slacks() const {
      return this->slacks;
   }

   const Collection<size_t>& AMPLModel::get_single_lower_bounded_variables() const {
      return this->single_lower_bounded_variables_collection;
   }

   const Collection<size_t>& AMPLModel::get_single_upper_bounded_variables() const {
      return this->single_upper_bounded_variables_collection;
   }

   const Collection<size_t>& AMPLModel::get_lower_bounded_variables() const {
      return this->lower_bounded_variables_collection;
   }

   const Collection<size_t>& AMPLModel::get_upper_bounded_variables() const {
      return this->upper_bounded_variables_collection;
   }

   bool are_all_zeros(const Vector<double>& multipliers) {
      return std::all_of(multipliers.begin(), multipliers.end(), [](double xj) {
         return xj == 0.;
      });
   }

   size_t AMPLModel::compute_hessian_number_nonzeros(double objective_multiplier, const Vector<double>& multipliers) const {
      // compute the sparsity
      const int objective_number = -1;
      const int upper_triangular = 1;
      const bool all_zeros_multipliers = are_all_zeros(multipliers);
      int number_nonzeros = (*(this->asl)->p.Sphset)(this->asl, nullptr, objective_number, (objective_multiplier != 0.),
            not all_zeros_multipliers, upper_triangular);
      return static_cast<size_t>(number_nonzeros);
   }

   void AMPLModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const {
      // register the vector of variables
      (*(this->asl)->p.Xknown)(this->asl, const_cast<double*>(x.data()), nullptr);

      // scale by the objective sign
      objective_multiplier *= this->objective_sign;

      // compute the number of nonzeros
      [[maybe_unused]] const size_t number_nonzeros = this->fixed_hessian_sparsity ? this->number_asl_hessian_nonzeros :
                                                      this->compute_hessian_number_nonzeros(objective_multiplier, multipliers);
      assert(hessian.capacity >= number_nonzeros);

      // evaluate the Hessian: store the matrix in a preallocated array this->asl_hessian
      const int objective_number = -1;
      if (this->fixed_hessian_sparsity) {
         (*(this->asl)->p.Sphes)(this->asl, nullptr, const_cast<double*>(this->asl_hessian.data()), objective_number, &objective_multiplier,
               const_cast<double*>(multipliers.data()));
      }
      else {
         double* objective_multiplier_pointer = (objective_multiplier != 0.) ? &objective_multiplier : nullptr;
         const bool all_zeros_multipliers = are_all_zeros(multipliers);
         (*(this->asl)->p.Sphes)(this->asl, nullptr, const_cast<double*>(this->asl_hessian.data()), objective_number, objective_multiplier_pointer,
               all_zeros_multipliers ? nullptr : const_cast<double*>(multipliers.data()));
      }

      // generate the sparsity pattern in the right sparse format
      const int* asl_column_start = this->asl->i.sputinfo_->hcolstarts;
      const int* asl_row_index = this->asl->i.sputinfo_->hrownos;
      // check that the column pointers are sorted in increasing order
      // TODO throw exception
      assert(in_increasing_order(asl_column_start, this->number_variables + 1) && "AMPLModel::evaluate_lagrangian_hessian: column starts are not ordered");

      // copy the nonzeros in the Hessian
      hessian.reset();
      for (size_t column_index: Range(this->number_variables)) {
         for (size_t k: Range(static_cast<size_t>(asl_column_start[column_index]), static_cast<size_t>(asl_column_start[column_index + 1]))) {
            const size_t row_index = static_cast<size_t>(asl_row_index[k]);
            const double entry = this->asl_hessian[k];
            hessian.insert(entry, row_index, column_index);
         }
         hessian.finalize_column(column_index);
      }
      // unregister the vector of variables
      this->asl->i.x_known = 0;
   }

   double AMPLModel::variable_lower_bound(size_t variable_index) const {
      return this->variable_lower_bounds[variable_index];
   }

   double AMPLModel::variable_upper_bound(size_t variable_index) const {
      return this->variable_upper_bounds[variable_index];
   }

   BoundType AMPLModel::get_variable_bound_type(size_t variable_index) const {
      return this->variable_status[variable_index];
   }

   double AMPLModel::constraint_lower_bound(size_t constraint_index) const {
      return this->constraint_lower_bounds[constraint_index];
   }

   double AMPLModel::constraint_upper_bound(size_t constraint_index) const {
      return this->constraint_upper_bounds[constraint_index];
   }

   FunctionType AMPLModel::get_constraint_type(size_t constraint_index) const {
      return this->constraint_type[constraint_index];
   }

   BoundType AMPLModel::get_constraint_bound_type(size_t constraint_index) const {
      return this->constraint_status[constraint_index];
   }

   // initial primal point
   void AMPLModel::initial_primal_point(Vector<double>& x) const {
      assert(x.size() >= this->number_variables);
      std::copy(this->asl->i.X0_, this->asl->i.X0_ + this->number_variables, x.begin());
   }

   // initial dual point
   void AMPLModel::initial_dual_point(Vector<double>& multipliers) const {
      assert(multipliers.size() >= this->number_constraints);
      std::copy(this->asl->i.pi0_, this->asl->i.pi0_ + this->number_constraints, multipliers.begin());
   }

   void AMPLModel::postprocess_solution(Iterate& /*iterate*/, TerminationStatus /*termination_status*/) const {
      // do nothing
   }

   void AMPLModel::generate_constraints() {
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->constraint_lower_bounds[constraint_index] = (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2*constraint_index] : -INF<double>;
         this->constraint_upper_bounds[constraint_index] = (this->asl->i.LUrhs_ != nullptr) ? this->asl->i.LUrhs_[2*constraint_index + 1] : INF<double>;
      }
      AMPLModel::determine_bounds_types(this->constraint_lower_bounds, this->constraint_upper_bounds, this->constraint_status);

      // partition equality and inequality constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (this->get_constraint_bound_type(constraint_index) == EQUAL_BOUNDS) {
            this->equality_constraints.push_back(constraint_index);
         }
         else {
            this->inequality_constraints.push_back(constraint_index);
         }
      }

      // AMPL orders the constraints based on the function type: nonlinear first, then linear
      const size_t number_nonlinear_constraints = static_cast<size_t>(this->asl->i.nlc_);
      for (size_t constraint_index: Range(number_nonlinear_constraints)) {
         this->constraint_type[constraint_index] = NONLINEAR;
      }
      for (size_t constraint_index: Range(number_nonlinear_constraints, this->number_constraints)) {
         this->constraint_type[constraint_index] = LINEAR;
         this->linear_constraints.push_back(constraint_index);
      }
   }

   void AMPLModel::determine_bounds_types(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds, std::vector<BoundType>& status) {
      assert(lower_bounds.size() == status.size());
      assert(upper_bounds.size() == status.size());
      // build the "status" vector as a mapping (map/transform operation) of the "bounds" vector
      for (size_t index: Range(lower_bounds.size())) {
         if (lower_bounds[index] == upper_bounds[index]) {
            status[index] = EQUAL_BOUNDS;
         }
         else if (is_finite(lower_bounds[index]) && is_finite(upper_bounds[index])) {
            status[index] = BOUNDED_BOTH_SIDES;
         }
         else if (is_finite(lower_bounds[index])) {
            status[index] = BOUNDED_LOWER;
         }
         else if (is_finite(upper_bounds[index])) {
            status[index] = BOUNDED_UPPER;
         }
         else {
            status[index] = UNBOUNDED;
         }
      }
   }
} // namespace
