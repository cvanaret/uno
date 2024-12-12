// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "OptimalityProblem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"

namespace uno {
   OptimalityProblem::OptimalityProblem(const Model& model): OptimizationProblem(model, model.number_variables, model.number_constraints) {
   }

   void OptimalityProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
      iterate.evaluate_objective_gradient(this->model);
      // TODO change this
      objective_gradient = iterate.evaluations.objective_gradient;
   }

   void OptimalityProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;
   }

   void OptimalityProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      iterate.evaluate_constraint_jacobian(this->model);
      // TODO change this
      constraint_jacobian = iterate.evaluations.constraint_jacobian;
   }

   void OptimalityProblem::evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      this->model.evaluate_lagrangian_hessian(x, this->get_objective_multiplier(), multipliers, hessian);
   }

   void OptimalityProblem::compute_hessian_vector_product(const Vector<double>& x, const Vector<double>& multipliers, Vector<double>& result) const {
      this->model.compute_hessian_vector_product(x, this->get_objective_multiplier(), multipliers, result);
   }

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void OptimalityProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const {
      lagrangian_gradient.objective_contribution.fill(0.);
      lagrangian_gradient.constraints_contribution.fill(0.);

      // objective gradient
      for (auto [variable_index, derivative]: iterate.evaluations.objective_gradient) {
         lagrangian_gradient.objective_contribution[variable_index] += derivative;
      }

      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         if (multipliers.constraints[constraint_index] != 0.) {
            for (auto [variable_index, derivative]: iterate.evaluations.constraint_jacobian[constraint_index]) {
               lagrangian_gradient.constraints_contribution[variable_index] -= multipliers.constraints[constraint_index] * derivative;
            }
         }
      }

      // bound constraints of original variables
      for (size_t variable_index: Range(this->number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (multipliers.lower_bounds[variable_index] +
                                                                          multipliers.upper_bounds[variable_index]);
      }
   }

   double OptimalityProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      // bound constraints
      const Range variables_range = Range(this->model.number_variables);
      const VectorExpression variable_complementarity{variables_range, [&](size_t variable_index) {
         if (0. < multipliers.lower_bounds[variable_index]) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->model.variable_lower_bound(variable_index)) - shift_value;
         }
         if (multipliers.upper_bounds[variable_index] < 0.) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->model.variable_upper_bound(variable_index)) - shift_value;
         }
         return 0.;
      }};

      // inequality constraints
      const VectorExpression constraint_complementarity{this->model.get_inequality_constraints(), [&](size_t constraint_index) {
         if (0. < multipliers.constraints[constraint_index]) { // lower bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->model.constraint_lower_bound(constraint_index)) -
                   shift_value;
         }
         else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
            return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->model.constraint_upper_bound(constraint_index)) -
                   shift_value;
         }
         return 0.;
      }};
      return norm(residual_norm, variable_complementarity, constraint_complementarity);
   }
} // namespace