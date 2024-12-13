// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimizationProblem.hpp"

namespace uno {
   OptimizationProblem::OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints):
         model(model), number_variables(number_variables), number_constraints(number_constraints) {
   }

   bool OptimizationProblem::is_constrained() const {
      return (0 < this->number_constraints);
   }

   bool OptimizationProblem::has_inequality_constraints() const {
      return (not this->model.get_inequality_constraints().empty());
   }

   bool OptimizationProblem::has_fixed_variables() const {
      return (not this->model.get_fixed_variables().empty());
   }

   size_t OptimizationProblem::get_number_original_variables() const {
      return this->model.number_variables;
   }

   double OptimizationProblem::stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
         Norm residual_norm) {
      // norm of the scaled Lagrangian gradient
      const auto scaled_lagrangian = objective_multiplier * lagrangian_gradient.objective_contribution + lagrangian_gradient.constraints_contribution;
      return norm(residual_norm, scaled_lagrangian);
   }
} // namespace