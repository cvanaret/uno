// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Iterate.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"

namespace uno {
   size_t Iterate::number_eval_objective = 0;
   size_t Iterate::number_eval_constraints = 0;
   size_t Iterate::number_eval_objective_gradient = 0;
   size_t Iterate::number_eval_jacobian = 0;

   Iterate::Iterate(size_t number_variables, size_t number_constraints) :
         number_variables(number_variables), number_constraints(number_constraints),
         primals(number_variables), multipliers(number_variables, number_constraints), feasibility_multipliers(number_variables, number_constraints),
         evaluations(number_variables, number_constraints), residuals(number_variables), feasibility_residuals(number_variables) {
   }

   void Iterate::evaluate_objective(const Model& model) {
      if (not this->is_objective_computed) {
         // evaluate the objective
         this->evaluations.objective = model.evaluate_objective(this->primals);
         Iterate::number_eval_objective++;
         if (not is_finite(this->evaluations.objective)) {
            throw FunctionEvaluationError();
         }
         this->is_objective_computed = true;
      }
   }

   void Iterate::evaluate_constraints(const Model& model) {
      if (not this->are_constraints_computed) {
         if (model.is_constrained()) {
            // evaluate the constraints
            model.evaluate_constraints(this->primals, this->evaluations.constraints);
            Iterate::number_eval_constraints++;
            // check finiteness
            if (std::any_of(this->evaluations.constraints.begin(), this->evaluations.constraints.end(), [](double constraint_j) {
               return not is_finite(constraint_j);
            })) {
               throw FunctionEvaluationError();
            }
         }
         this->are_constraints_computed = true;
      }
   }

   void Iterate::evaluate_objective_gradient(const Model& model) {
      if (not this->is_objective_gradient_computed) {
         this->evaluations.objective_gradient.clear();
         // evaluate the objective gradient
         model.evaluate_objective_gradient(this->primals, this->evaluations.objective_gradient);
         this->is_objective_gradient_computed = true;
         Iterate::number_eval_objective_gradient++;
      }
   }

   void Iterate::evaluate_constraint_jacobian(const Model& model) {
      if (not this->is_constraint_jacobian_computed) {
         this->evaluations.constraint_jacobian.clear();
         if (model.is_constrained()) {
            model.evaluate_constraint_jacobian(this->primals, this->evaluations.constraint_jacobian);
            Iterate::number_eval_jacobian++;
         }
         this->is_constraint_jacobian_computed = true;
      }
   }

   void Iterate::set_number_variables(size_t new_number_variables) {
      this->number_variables = new_number_variables;
      this->primals.resize(new_number_variables);
      this->multipliers.lower_bounds.resize(new_number_variables);
      this->multipliers.upper_bounds.resize(new_number_variables);
      this->evaluations.objective_gradient.reserve(new_number_variables);
      this->residuals.lagrangian_gradient.resize(new_number_variables);
   }

   std::ostream& operator<<(std::ostream& stream, const Iterate& iterate) {
      stream << "Primal variables: " << iterate.primals << '\n';
      stream << "            ┌ Constraint: " << iterate.multipliers.constraints << '\n';
      stream << "Multipliers │ Lower bound: " << iterate.multipliers.lower_bounds << '\n';
      stream << "            └ Upper bound: " << iterate.multipliers.upper_bounds << '\n';
      stream << "                        ┌ Constraint: " << iterate.feasibility_multipliers.constraints << '\n';
      stream << "Feasibility multipliers │ Lower bound: " << iterate.feasibility_multipliers.lower_bounds << '\n';
      stream << "                        └ Upper bound: " << iterate.feasibility_multipliers.upper_bounds << '\n';
      stream << "Objective value: " << iterate.evaluations.objective << '\n';
      stream << "Primal feasibility: " << iterate.primal_feasibility << '\n';

      stream << "          ┌ Stationarity: " << iterate.residuals.stationarity << '\n';
      stream << "Residuals │ Complementarity: " << iterate.residuals.complementarity << '\n';
      stream << "          └ Lagrangian gradient: " << iterate.residuals.lagrangian_gradient;
      stream << "Feasibility residuals ┌ Stationarity: " << iterate.feasibility_residuals.stationarity << '\n';
      stream << "                      │ Complementarity: " << iterate.feasibility_residuals.complementarity << '\n';
      stream << "                      └ Lagrangian gradient: " << iterate.feasibility_residuals.lagrangian_gradient;

      stream << "                  ┌ Infeasibility: " << iterate.progress.infeasibility << '\n';
      stream << "Progress measures │ Optimality: " << iterate.progress.objective(1.) << '\n';
      stream << "                  └ Auxiliary terms: " << iterate.progress.auxiliary << '\n';


      return stream;
   }
} // namespace