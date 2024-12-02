// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ScaledModel.hpp"
#include "Model.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace uno {
   ScaledModel::ScaledModel(std::unique_ptr<Model> original_model, Iterate& initial_iterate, const Options& options):
         Model(original_model->name + " -> scaled", original_model->number_variables, original_model->number_constraints,
               original_model->objective_sign),
         model(std::move(original_model)),
         scaling(this->model->number_constraints, options.get_double("function_scaling_threshold")),
         scaled_multipliers(this->number_constraints) {
      if (options.get_bool("scale_functions")) {
         // evaluate the gradients at the current point
         initial_iterate.evaluate_objective_gradient(*this->model);
         initial_iterate.evaluate_constraint_jacobian(*this->model);
         this->scaling.compute(initial_iterate.evaluations.objective_gradient, initial_iterate.evaluations.constraint_jacobian);
         // scale the gradients
         scale(initial_iterate.evaluations.objective_gradient, this->scaling.get_objective_scaling());
         for (size_t constraint_index: Range(this->model->number_constraints)) {
            scale(initial_iterate.evaluations.constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
         }
         // since the definition of the constraints changed, reset the evaluation flags
         initial_iterate.is_objective_gradient_computed = false;
         initial_iterate.is_constraint_jacobian_computed = false;
      }
      // check the scaling factors
      assert(0 < this->scaling.get_objective_scaling() && "Objective scaling failed.");
      for ([[maybe_unused]] size_t constraint_index: Range(this->number_constraints)) {
         assert(0 < this->scaling.get_constraint_scaling(constraint_index) && "Constraint scaling failed.");
      }
   }

   double ScaledModel::evaluate_objective(const Vector<double>& x) const {
      const double objective = this->model->evaluate_objective(x);
      return this->scaling.get_objective_scaling()*objective;
   }

   void ScaledModel::evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const {
      this->model->evaluate_objective_gradient(x, gradient);
      scale(gradient, this->scaling.get_objective_scaling());
   }

   void ScaledModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      this->model->evaluate_constraints(x, constraints);
      for (size_t constraint_index: Range(this->number_constraints)) {
         constraints[constraint_index] *= this->scaling.get_constraint_scaling(constraint_index);
      }
   }

   void ScaledModel::evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
      this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
      scale(gradient, this->scaling.get_constraint_scaling(constraint_index));
   }

   void ScaledModel::evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
      this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
      for (size_t constraint_index: Range(this->number_constraints)) {
         scale(constraint_jacobian[constraint_index], this->scaling.get_constraint_scaling(constraint_index));
      }
   }

   void ScaledModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      // scale the objective and constraint multipliers
      const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->scaled_multipliers[constraint_index] = this->scaling.get_constraint_scaling(constraint_index) * multipliers[constraint_index];
      }
      this->model->evaluate_lagrangian_hessian(x, scaled_objective_multiplier, this->scaled_multipliers, hessian);
   }

   void ScaledModel::compute_hessian_vector_product(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         Vector<double>& result) const {
      // scale the objective and constraint multipliers
      const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
      for (size_t constraint_index: Range(this->number_constraints)) {
         this->scaled_multipliers[constraint_index] = this->scaling.get_constraint_scaling(constraint_index) * multipliers[constraint_index];
      }
      this->model->compute_hessian_vector_product(x, scaled_objective_multiplier, this->scaled_multipliers, result);
   }

   double ScaledModel::variable_lower_bound(size_t variable_index) const {
      return this->model->variable_lower_bound(variable_index);
   }

   double ScaledModel::variable_upper_bound(size_t variable_index) const {
      return this->model->variable_upper_bound(variable_index);
   }

   BoundType ScaledModel::get_variable_bound_type(size_t variable_index) const {
      return this->model->get_variable_bound_type(variable_index);
   }

   const Collection<size_t>& ScaledModel::get_lower_bounded_variables() const {
      return this->model->get_lower_bounded_variables();
   }

   const Collection<size_t>& ScaledModel::get_upper_bounded_variables() const {
      return this->model->get_upper_bounded_variables();
   }

   const SparseVector<size_t>& ScaledModel::get_slacks() const {
      return this->model->get_slacks();
   }

   const Collection<size_t>& ScaledModel::get_single_lower_bounded_variables() const {
      return this->model->get_single_lower_bounded_variables();
   }

   const Collection<size_t>& ScaledModel::get_single_upper_bounded_variables() const {
      return this->model->get_single_upper_bounded_variables();
   }

   const Vector<size_t>& ScaledModel::get_fixed_variables() const {
      return this->model->get_fixed_variables();
   }

   double ScaledModel::constraint_lower_bound(size_t constraint_index) const {
      return this->scaling.get_constraint_scaling(constraint_index) * this->model->constraint_lower_bound(constraint_index);
   }

   double ScaledModel::constraint_upper_bound(size_t constraint_index) const {
      return this->scaling.get_constraint_scaling(constraint_index) * this->model->constraint_upper_bound(constraint_index);
   }

   FunctionType ScaledModel::get_constraint_type(size_t constraint_index) const {
      return this->model->get_constraint_type(constraint_index);
   }

   BoundType ScaledModel::get_constraint_bound_type(size_t constraint_index) const {
      return this->model->get_constraint_bound_type(constraint_index);
   }

   const Collection<size_t>& ScaledModel::get_equality_constraints() const {
      return this->model->get_equality_constraints();
   }

   const Collection<size_t>& ScaledModel::get_inequality_constraints() const {
      return this->model->get_inequality_constraints();
   }

   const Collection<size_t>& ScaledModel::get_linear_constraints() const {
      return this->model->get_linear_constraints();
   }

   void ScaledModel::initial_primal_point(Vector<double>& x) const {
      this->model->initial_primal_point(x);
   }
   
   void ScaledModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model->initial_dual_point(multipliers);
   }

   void ScaledModel::postprocess_solution(Iterate& iterate, IterateStatus termination_status) const {
      // unscale the objective value
      if (iterate.is_objective_computed) {
         iterate.evaluations.objective /= this->scaling.get_objective_scaling();
      }

      // unscale the constraint multipliers
      for (size_t constraint_index: Range(iterate.number_constraints)) {
         iterate.multipliers.constraints[constraint_index] *= this->scaling.get_constraint_scaling(constraint_index) / this->scaling.get_objective_scaling();
      }

      // unscale the bound multipliers
      for (size_t variable_index: Range(iterate.number_variables)) {
         iterate.multipliers.lower_bounds[variable_index] /= this->scaling.get_objective_scaling();
         iterate.multipliers.upper_bounds[variable_index] /= this->scaling.get_objective_scaling();
      }
      this->model->postprocess_solution(iterate, termination_status);
   }

   size_t ScaledModel::number_objective_gradient_nonzeros() const {
      return this->model->number_objective_gradient_nonzeros();
   }

   size_t ScaledModel::number_jacobian_nonzeros() const {
      return this->model->number_jacobian_nonzeros();
   }

   size_t ScaledModel::number_hessian_nonzeros() const {
      return this->model->number_hessian_nonzeros();
   }
} // namespace

