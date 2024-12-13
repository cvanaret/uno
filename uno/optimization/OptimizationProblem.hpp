// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include "linear_algebra/Norm.hpp"
#include "model/Model.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   class Iterate;
   class Multipliers;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   class OptimizationProblem {
   public:
      OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints);
      virtual ~OptimizationProblem() = default;

      const Model& model;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */

      [[nodiscard]] bool is_constrained() const;
      [[nodiscard]] bool has_inequality_constraints() const;
      [[nodiscard]] bool has_fixed_variables() const;

      // function evaluations
      [[nodiscard]] virtual double get_objective_multiplier() const = 0;
      virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const = 0;
      virtual void evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const = 0;
      virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const = 0;
      virtual void evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers, SymmetricMatrix<size_t, double>& hessian) const = 0;

      [[nodiscard]] size_t get_number_original_variables() const;
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_lower_bounded_variables() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_upper_bounded_variables() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_single_lower_bounded_variables() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_single_upper_bounded_variables() const = 0;

      [[nodiscard]] virtual size_t number_objective_gradient_nonzeros() const = 0;
      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const = 0;
      [[nodiscard]] virtual size_t number_hessian_nonzeros() const = 0;

      [[nodiscard]] static double stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
            Norm residual_norm);
      virtual void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate, const Multipliers& multipliers) const = 0;
      [[nodiscard]] virtual double complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
            const Multipliers& multipliers, double shift_value, Norm residual_norm) const = 0;
   };
} // namespace

#endif // UNO_OPTIMIZATIONPROBLEM_H