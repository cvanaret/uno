// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include "optimization/OptimizationProblem.hpp"
#include "ElasticVariables.hpp"
#include "symbolic/Concatenation.hpp"

namespace uno {
   class l1RelaxedProblem: public OptimizationProblem {
   public:
      l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient, double proximal_coefficient,
            double const* proximal_center);

      [[nodiscard]] double get_objective_multiplier() const override;
      void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
      void evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const override;
      void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers, SymmetricMatrix<size_t, double>& hessian) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

      void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate, const Multipliers& multipliers) const override;
      [[nodiscard]] double complementarity_error(const Vector<double>& primals, const Vector<double>& constraints, const Multipliers& multipliers,
            double shift_value, Norm residual_norm) const override;

      // parameterization
      void set_objective_multiplier(double new_objective_multiplier);

      void set_proximal_multiplier(double new_proximal_coefficient);
      void set_proximal_center(double const* new_proximal_center);
      void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

   protected:
      double objective_multiplier;
      const double constraint_violation_coefficient;
      double proximal_coefficient;
      double const* proximal_center;
      ElasticVariables elastic_variables;
      const Concatenation<const Collection<size_t>&, ForwardRange> lower_bounded_variables; // model variables + elastic variables
      const Concatenation<const Collection<size_t>&, ForwardRange> single_lower_bounded_variables; // model variables + elastic variables

      // delegating constructor
      l1RelaxedProblem(const Model& model, ElasticVariables&& elastic_variables, double objective_multiplier, double constraint_violation_coefficient,
            double proximal_coefficient, double const* proximal_center);
   };
} // namespace

#endif // UNO_L1RELAXEDPROBLEM_H
