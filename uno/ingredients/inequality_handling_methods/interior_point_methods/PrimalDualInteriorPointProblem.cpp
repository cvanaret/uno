// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   PrimalDualInteriorPointProblem::PrimalDualInteriorPointProblem(const OptimizationProblem& problem, const Multipliers& current_multipliers,
         double barrier_parameter):
         OptimizationProblem(problem.model, problem.number_variables, problem.number_constraints),
         problem(problem), current_multipliers(current_multipliers), barrier_parameter(barrier_parameter) { }

   double PrimalDualInteriorPointProblem::get_objective_multiplier() const {
      return this->problem.get_objective_multiplier();
   }

   void PrimalDualInteriorPointProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
      this->problem.evaluate_objective_gradient(iterate, objective_gradient);

      // barrier terms
      for (size_t variable_index: Range(this->problem.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(problem.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - problem.variable_lower_bound(variable_index));
            // damping
            if (not is_finite(problem.variable_upper_bound(variable_index))) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(problem.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - problem.variable_upper_bound(variable_index));
            // damping
            if (not is_finite(problem.variable_lower_bound(variable_index))) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         objective_gradient.insert(variable_index, barrier_term);
      }
   }

   void PrimalDualInteriorPointProblem::evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const {
      problem.evaluate_constraints(iterate, constraints);
   }

   void PrimalDualInteriorPointProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
      problem.evaluate_constraint_jacobian(iterate, constraint_jacobian);
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      // original Lagrangian Hessian
      this->problem.evaluate_lagrangian_hessian(x, multipliers, hessian);

      // barrier terms
      for (size_t variable_index: Range(this->problem.number_variables)) {
         double diagonal_barrier_term = 0.;
         if (is_finite(this->problem.variable_lower_bound(variable_index))) { // lower bounded
            const double distance_to_bound = x[variable_index] - problem.variable_lower_bound(variable_index);
            diagonal_barrier_term += this->current_multipliers.lower_bounds[variable_index] / distance_to_bound;
         }
         if (is_finite(this->problem.variable_upper_bound(variable_index))) { // upper bounded
            const double distance_to_bound = x[variable_index] - problem.variable_upper_bound(variable_index);
            diagonal_barrier_term += this->current_multipliers.upper_bounds[variable_index] / distance_to_bound;
         }
         hessian.insert(diagonal_barrier_term, variable_index, variable_index);
      }
   }

   double PrimalDualInteriorPointProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalDualInteriorPointProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   double PrimalDualInteriorPointProblem::constraint_lower_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   double PrimalDualInteriorPointProblem::constraint_upper_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_lower_bounded_variables() const {
      return this->problem.get_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_upper_bounded_variables() const {
      return this->problem.get_upper_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_lower_bounded_variables() const {
      return this->problem.get_single_lower_bounded_variables();
   }

   const Collection<size_t>& PrimalDualInteriorPointProblem::get_single_upper_bounded_variables() const {
      return this->problem.get_single_upper_bounded_variables();
   }

   size_t PrimalDualInteriorPointProblem::number_objective_gradient_nonzeros() const {
      size_t number_nonzeros = this->problem.number_objective_gradient_nonzeros();
      // barrier contribution
      for (size_t variable_index: Range(this->problem.number_variables)) {
         if (is_finite(this->problem.variable_lower_bound(variable_index)) || is_finite(this->problem.variable_upper_bound(variable_index))) {
            number_nonzeros++;
         }
      }
      return number_nonzeros;
   }

   size_t PrimalDualInteriorPointProblem::number_jacobian_nonzeros() const {
      return this->problem.number_jacobian_nonzeros();
   }

   size_t PrimalDualInteriorPointProblem::number_hessian_nonzeros() const {
      size_t number_zeros = this->problem.number_hessian_nonzeros();
      // barrier contribution
      for (size_t variable_index: Range(this->problem.number_variables)) {
         if (is_finite(this->problem.variable_lower_bound(variable_index)) || is_finite(this->problem.variable_upper_bound(variable_index))) {
            number_zeros++;
         }
      }
      return number_zeros;
   }

   void PrimalDualInteriorPointProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const {
      this->problem.evaluate_lagrangian_gradient(lagrangian_gradient, iterate, multipliers);

      // barrier terms
      for (size_t variable_index: Range(this->problem.number_variables)) {
         double barrier_term = 0.;
         if (is_finite(problem.variable_lower_bound(variable_index))) { // lower bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - problem.variable_lower_bound(variable_index));
            // damping
            if (not is_finite(problem.variable_upper_bound(variable_index))) {
               barrier_term += this->damping_factor * this->barrier_parameter;
            }
         }
         if (is_finite(problem.variable_upper_bound(variable_index))) { // upper bounded
            barrier_term += -this->barrier_parameter/(iterate.primals[variable_index] - problem.variable_upper_bound(variable_index));
            // damping
            if (not is_finite(problem.variable_lower_bound(variable_index))) {
               barrier_term -= this->damping_factor * this->barrier_parameter;
            }
         }
         // the objective contribution of the Lagrangian gradient may be scaled. Barrier terms go into the constraint contribution
         lagrangian_gradient.constraints_contribution[variable_index] += barrier_term;
      }
   }

   double PrimalDualInteriorPointProblem::complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->problem.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }
} // namespace