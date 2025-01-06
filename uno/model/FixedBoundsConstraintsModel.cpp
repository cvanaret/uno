// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedBoundsConstraintsModel.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   FixedBoundsConstraintsModel::FixedBoundsConstraintsModel(std::unique_ptr<Model> original_model, const Options& /*options*/) :
         Model(original_model->name + " -> no fixed bounds", original_model->number_variables,
               // move the fixed variables to the set of general constraints
               original_model->number_constraints + original_model->get_fixed_variables().size(),
               original_model->objective_sign),
         model(std::move(original_model)),
         fixed_variables(), // empty vector
         lower_bounded_variables_collection(this->lower_bounded_variables),
         upper_bounded_variables_collection(this->upper_bounded_variables),
         equality_constraints(concatenate(this->model->get_equality_constraints(), Range(this->model->number_constraints, this->number_constraints))),
         linear_constraints(concatenate(this->model->get_linear_constraints(), Range(this->model->number_constraints, this->number_constraints))) {
      this->lower_bounded_variables.reserve(this->model->get_lower_bounded_variables().size());
      this->upper_bounded_variables.reserve(this->model->get_upper_bounded_variables().size());

      for (size_t variable_index: Range(this->model->number_variables)) {
         if (is_finite(this->variable_lower_bound(variable_index))) {
            this->lower_bounded_variables.emplace_back(variable_index);
         }
         if (is_finite(this->variable_upper_bound(variable_index))) {
            this->upper_bounded_variables.emplace_back(variable_index);
         }
      }
   }

   double FixedBoundsConstraintsModel::evaluate_objective(const Vector<double>& x) const {
      return this->model->evaluate_objective(x);
   }

   void FixedBoundsConstraintsModel::evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const {
      this->model->evaluate_objective_gradient(x, gradient);
   }

   void FixedBoundsConstraintsModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->model->evaluate_constraints(x, constraints);
      // add the fixed variables
      size_t current_constraint = this->model->number_constraints;
      for (size_t fixed_variable_index: this->model->get_fixed_variables()) {
         constraints[current_constraint] = x[fixed_variable_index];
         current_constraint++;
      }
   }

   void FixedBoundsConstraintsModel::evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
      if (constraint_index < this->model->number_constraints) {
         // original constraint
         this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
      }
      else {
         // fixed variable
         const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
         gradient.insert(variable_index, 1.);
      }
   }

   void FixedBoundsConstraintsModel::evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
      this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
      // add the fixed variables
      size_t current_constraint = this->model->number_constraints;
      for (size_t fixed_variable_index: this->model->get_fixed_variables()) {
         constraint_jacobian[current_constraint].insert(fixed_variable_index, 1.);
         current_constraint++;
      }
   }

   void FixedBoundsConstraintsModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   }

   double FixedBoundsConstraintsModel::variable_lower_bound(size_t variable_index) const {
      if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
      // remove bounds of fixed variables
         return -INF<double>;
      }
      return this->model->variable_lower_bound(variable_index);
   }

   double FixedBoundsConstraintsModel::variable_upper_bound(size_t variable_index) const {
      if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
      // remove bounds of fixed variables
         return INF<double>;
      }
      return this->model->variable_upper_bound(variable_index);
   }

   BoundType FixedBoundsConstraintsModel::get_variable_bound_type(size_t variable_index) const {
      if (this->model->variable_lower_bound(variable_index) == this->model->variable_upper_bound(variable_index)) {
      // fixed variable: remove the bounds
         return BoundType::UNBOUNDED;
      }
      else {
         return this->model->get_variable_bound_type(variable_index);
      }
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_lower_bounded_variables() const {
      return this->lower_bounded_variables_collection;
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_upper_bounded_variables() const {
      return this->upper_bounded_variables_collection;
   }

   const SparseVector<size_t>& FixedBoundsConstraintsModel::get_slacks() const {
      return this->model->get_slacks();
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_single_lower_bounded_variables() const {
      return this->model->get_single_lower_bounded_variables();
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_single_upper_bounded_variables() const {
      return this->model->get_single_upper_bounded_variables();
   }

   const Vector<size_t>& FixedBoundsConstraintsModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double FixedBoundsConstraintsModel::constraint_lower_bound(size_t constraint_index) const {
      if (constraint_index < this->model->number_constraints) {
// original constraint
         return this->model->constraint_lower_bound(constraint_index);
      }
      else {
// fixed variable
         const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
         return this->model->variable_lower_bound(variable_index);
      }
   }

   double FixedBoundsConstraintsModel::constraint_upper_bound(size_t constraint_index) const {
      if (constraint_index < this->model->number_constraints) {
// original constraint
         return this->model->constraint_upper_bound(constraint_index);
      }
      else {
// fixed variable
         const size_t variable_index = this->model->get_fixed_variables()[constraint_index - this->model->number_constraints];
         return this->model->variable_lower_bound(variable_index);
      }
   }

   FunctionType FixedBoundsConstraintsModel::get_constraint_type(size_t constraint_index) const {
      if (constraint_index < this->model->number_constraints) {
// original constraint
         return this->model->get_constraint_type(constraint_index);
      }
      else {
// fixed variables
         return FunctionType::LINEAR;
      }
   }

   BoundType FixedBoundsConstraintsModel::get_constraint_bound_type(size_t constraint_index) const {
      if (constraint_index < this->model->number_constraints) {
// original constraint
         return this->model->get_constraint_bound_type(constraint_index);
      }
      else {
// fixed variables
         return BoundType::EQUAL_BOUNDS;
      }
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_equality_constraints() const {
      return this->equality_constraints;
   }
   const Collection<size_t>& FixedBoundsConstraintsModel::get_inequality_constraints() const {
      return this->model->get_inequality_constraints();
   }
   const Collection<size_t>& FixedBoundsConstraintsModel::get_linear_constraints() const {
      return this->linear_constraints;
   }

   void FixedBoundsConstraintsModel::initial_primal_point(Vector<double>& x) const {
      this->model->initial_primal_point(x);
// set the fixed variables
      for (size_t variable_index: this->model->get_fixed_variables()) {
         x[variable_index] = this->model->variable_lower_bound(variable_index);
      }
   }

   void FixedBoundsConstraintsModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model->initial_dual_point(multipliers);
   }

   void FixedBoundsConstraintsModel::postprocess_solution(Iterate& iterate, IterateStatus termination_status) const {
// move the multipliers back from the general constraints to the bound constraints
      size_t current_constraint = this->model->number_constraints;
      for (size_t variable_index: this->model->get_fixed_variables()) {
         const double constraint_multiplier = iterate.multipliers.constraints[current_constraint];
         if (0. < constraint_multiplier) {
            iterate.multipliers.lower_bounds[variable_index] = constraint_multiplier;
         }
         else {
            iterate.multipliers.upper_bounds[variable_index] = constraint_multiplier;
         }
         current_constraint++;
      }
      this->model->postprocess_solution(iterate, termination_status);
   }

   size_t FixedBoundsConstraintsModel::number_objective_gradient_nonzeros() const {
      return this->model->number_objective_gradient_nonzeros();
   }

   size_t FixedBoundsConstraintsModel::number_jacobian_nonzeros() const {
      return this->model->number_jacobian_nonzeros() + this->model->get_fixed_variables().size();
   }

   size_t FixedBoundsConstraintsModel::number_hessian_nonzeros() const {
      return this->model->number_hessian_nonzeros();
   }
} // namespace