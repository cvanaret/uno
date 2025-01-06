// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "l1RelaxedProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorExpression.hpp"
#include "symbolic/Concatenation.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient,
         double proximal_coefficient, double const* proximal_center):
   // call delegating constructor
         l1RelaxedProblem(model, ElasticVariables::generate(model), objective_multiplier, constraint_violation_coefficient, proximal_coefficient,
               proximal_center) {
   }

   // private delegating constructor
   l1RelaxedProblem::l1RelaxedProblem(const Model& model, ElasticVariables&& elastic_variables, double objective_multiplier,
         double constraint_violation_coefficient, double proximal_coefficient, double const* proximal_center):
         OptimizationProblem(model, model.number_variables + elastic_variables.size(), model.number_constraints),
         objective_multiplier(objective_multiplier),
         constraint_violation_coefficient(constraint_violation_coefficient),
         proximal_coefficient(proximal_coefficient),
         proximal_center(proximal_center),
         elastic_variables(std::forward<ElasticVariables>(elastic_variables)),
         // lower bounded variables are the model variables + the elastic variables
         lower_bounded_variables(concatenate(this->model.get_lower_bounded_variables(), Range(model.number_variables,
               model.number_variables + this->elastic_variables.size()))),
         single_lower_bounded_variables(concatenate(this->model.get_single_lower_bounded_variables(),
               Range(model.number_variables, model.number_variables + this->elastic_variables.size()))) {
   }

   double l1RelaxedProblem::get_objective_multiplier() const {
      return this->objective_multiplier;
   }

   void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
      // scale nabla f(x) by rho
      if (this->objective_multiplier != 0.) {
         iterate.evaluate_objective_gradient(this->model);
         // TODO change this
         objective_gradient = iterate.evaluations.objective_gradient;
         scale(objective_gradient, this->objective_multiplier);
      }
      else {
         objective_gradient.clear();
      }

      // constraint violation (through elastic variables) contribution
      for (const auto [_, elastic_index]: this->elastic_variables.positive) {
         objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
      }
      for (const auto [_, elastic_index]: this->elastic_variables.negative) {
         objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling * (iterate.primals[variable_index] - this->proximal_center[variable_index]);
            objective_gradient.insert(variable_index, proximal_term);
         }
      }
   }

   void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const {
      iterate.evaluate_constraints(this->model);
      constraints = iterate.evaluations.constraints;
      // add the contribution of the elastics
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         constraints[constraint_index] -= iterate.primals[elastic_index];
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         constraints[constraint_index] += iterate.primals[elastic_index];
      }
   }

   void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      iterate.evaluate_constraint_jacobian(this->model);
      // TODO change this
      constraint_jacobian = iterate.evaluations.constraint_jacobian;
      // add the contribution of the elastics
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         constraint_jacobian[constraint_index].insert(elastic_index, -1.);
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         constraint_jacobian[constraint_index].insert(elastic_index, 1.);
      }
   }

   void l1RelaxedProblem::evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      this->model.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            hessian.insert(proximal_term, variable_index, variable_index);
         }
      }

      // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
      for (size_t constraint_index: Range(this->model.number_variables, this->number_variables)) {
         hessian.finalize_column(constraint_index);
      }
   }

   // Lagrangian gradient split in two parts: objective contribution and constraints' contribution
   void l1RelaxedProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
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
      for (size_t variable_index: Range(this->model.number_variables)) {
         lagrangian_gradient.constraints_contribution[variable_index] -= (multipliers.lower_bounds[variable_index] +
                                                                          multipliers.upper_bounds[variable_index]);
      }

      // elastic variables
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient +
            multipliers.constraints[constraint_index] - multipliers.lower_bounds[elastic_index];
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         lagrangian_gradient.constraints_contribution[elastic_index] += this->constraint_violation_coefficient -
            multipliers.constraints[constraint_index] - multipliers.lower_bounds[elastic_index];
      }

      // proximal contribution
      if (this->proximal_center != nullptr && this->proximal_coefficient != 0.) {
         for (size_t variable_index: Range(this->model.number_variables)) {
            const double scaling = std::min(1., 1./std::abs(this->proximal_center[variable_index]));
            const double proximal_term = this->proximal_coefficient * scaling * scaling;
            lagrangian_gradient.constraints_contribution[variable_index] += proximal_term * (iterate.primals[variable_index] -
                                                                                             this->proximal_center[variable_index]);
         }
      }
   }

   // complementary slackness error: expression for violated constraints depends on the definition of the relaxed problem
   double l1RelaxedProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& constraints, const Multipliers& multipliers,
         double shift_value, Norm residual_norm) const {
      // bound constraints
      const Range variables_range = Range(this->number_variables);
      const VectorExpression bounds_complementarity{variables_range, [&](size_t variable_index) {
         if (0. < multipliers.lower_bounds[variable_index]) {
            return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->variable_lower_bound(variable_index)) - shift_value;
         }
         if (multipliers.upper_bounds[variable_index] < 0.) {
            return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->variable_upper_bound(variable_index)) - shift_value;
         }
         return 0.;
      }};

      // general constraints
      // TODO use the values of the relaxed constraints
      const Range constraints_range = Range(this->number_constraints);
      const VectorExpression constraints_complementarity{constraints_range, [&](size_t constraint_index) {
         if (this->model.get_constraint_bound_type(constraint_index) != EQUAL_BOUNDS) {
            if (0. < multipliers.constraints[constraint_index]) { // lower bound
               return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_lower_bound(constraint_index)) - shift_value;
            }
            else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
               return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_upper_bound(constraint_index)) - shift_value;
            }
         }
         return 0.;
      }};
      return norm(residual_norm, bounds_complementarity);
   }

   double l1RelaxedProblem::variable_lower_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_lower_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return 0.;
      }
   }

   double l1RelaxedProblem::variable_upper_bound(size_t variable_index) const {
      if (variable_index < this->model.number_variables) { // model variable
         return this->model.variable_upper_bound(variable_index);
      }
      else { // elastic variable in [0, +inf[
         return INF<double>;
      }
   }

   double l1RelaxedProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->model.constraint_lower_bound(constraint_index);
   }

   double l1RelaxedProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->model.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& l1RelaxedProblem::get_lower_bounded_variables() const {
      return this->lower_bounded_variables;
   }

   const Collection<size_t>& l1RelaxedProblem::get_upper_bounded_variables() const {
      // same set as the model
      return this->model.get_upper_bounded_variables();
   }

   const Collection<size_t>& l1RelaxedProblem::get_single_lower_bounded_variables() const {
      return this->single_lower_bounded_variables;
   }

   const Collection<size_t>& l1RelaxedProblem::get_single_upper_bounded_variables() const {
      // same set as the model
      return this->model.get_single_upper_bounded_variables();
   }

   size_t l1RelaxedProblem::number_objective_gradient_nonzeros() const {
      // elastic contribution
      size_t number_nonzeros = this->elastic_variables.size();

      // objective contribution
      if (this->objective_multiplier != 0.) {
         number_nonzeros += this->model.number_objective_gradient_nonzeros();
      }
      return number_nonzeros;
   }

   size_t l1RelaxedProblem::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros() + this->elastic_variables.size();
   }

   size_t l1RelaxedProblem::number_hessian_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   void l1RelaxedProblem::set_objective_multiplier(double new_objective_multiplier) {
      assert(0. <= new_objective_multiplier && "The objective multiplier should be non-negative");
      this->objective_multiplier = new_objective_multiplier;
   }

   void l1RelaxedProblem::set_proximal_multiplier(double new_proximal_coefficient) {
      this->proximal_coefficient = new_proximal_coefficient;
   }

   void l1RelaxedProblem::set_proximal_center(double const* new_proximal_center) {
      this->proximal_center = new_proximal_center;
   }

   void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>&
   elastic_setting_function) const {
      iterate.set_number_variables(this->number_variables);
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
         elastic_setting_function(iterate, constraint_index, elastic_index, -1.);
      }
      for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
         elastic_setting_function(iterate, constraint_index, elastic_index, 1.);
      }
   }
} // namespace