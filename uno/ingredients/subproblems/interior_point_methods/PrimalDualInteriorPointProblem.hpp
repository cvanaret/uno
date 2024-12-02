// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTPROBLEM_H
#define UNO_PRIMALDUALINTERIORPOINTPROBLEM_H

#include "reformulation/OptimizationProblem.hpp"

namespace uno {
   class PrimalDualInteriorPointProblem : public OptimizationProblem {
   public:
      PrimalDualInteriorPointProblem(const OptimizationProblem& problem, const Multipliers& current_multipliers, double barrier_parameter,
            double damping_factor);

      // function evaluations
      [[nodiscard]] double get_objective_multiplier() const override;
      void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
      void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
      void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const override;
      void compute_hessian_vector_product(const Vector<double>& x, const Vector<double>& multipliers, Vector<double>& result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

      void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
            const Multipliers& multipliers) const override;
      [[nodiscard]] double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
            const Multipliers& multipliers, double shift_value, Norm residual_norm) const override;

   protected:
      const OptimizationProblem& problem;
      const Multipliers& current_multipliers;
      const double barrier_parameter;
      const double damping_factor;
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTPROBLEM_H
