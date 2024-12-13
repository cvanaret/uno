// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODEL_H
#define UNO_MODEL_H

#include <string>
#include "linear_algebra/Norm.hpp"
#include "optimization/IterateStatus.hpp"
#include "symbolic/VectorExpression.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   class Iterate;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   enum FunctionType {LINEAR, NONLINEAR};
   enum BoundType {EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED};

   class Model {
   public:
      Model(std::string name, size_t number_variables, size_t number_constraints, double objective_sign);
      virtual ~Model() = default;

      const std::string name;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */
      const double objective_sign; /*!< Sign of the objective function (1: minimization, -1: maximization) */

      [[nodiscard]] virtual double evaluate_objective(const Vector<double>& x) const = 0;
      virtual void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const = 0;
      virtual void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const = 0;
      virtual void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const = 0;
      virtual void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const = 0;
      virtual void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const = 0;

      // purely virtual functions
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const = 0;
      [[nodiscard]] virtual BoundType get_variable_bound_type(size_t variable_index) const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_lower_bounded_variables() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_upper_bounded_variables() const = 0;
      [[nodiscard]] virtual const SparseVector<size_t>& get_slacks() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_single_lower_bounded_variables() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_single_upper_bounded_variables() const = 0;
      [[nodiscard]] virtual const Vector<size_t>& get_fixed_variables() const = 0;

      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const = 0;
      [[nodiscard]] virtual FunctionType get_constraint_type(size_t constraint_index) const = 0;
      [[nodiscard]] virtual BoundType get_constraint_bound_type(size_t constraint_index) const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_equality_constraints() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const = 0;
      [[nodiscard]] virtual const Collection<size_t>& get_linear_constraints() const = 0;

      virtual void initial_primal_point(Vector<double>& x) const = 0;
      virtual void initial_dual_point(Vector<double>& multipliers) const = 0;
      virtual void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const = 0;

      [[nodiscard]] virtual size_t number_objective_gradient_nonzeros() const = 0;
      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const = 0;
      [[nodiscard]] virtual size_t number_hessian_nonzeros() const = 0;

      // auxiliary functions
      void project_onto_variable_bounds(Vector<double>& x) const;
      [[nodiscard]] bool is_constrained() const;

      // constraint violation
      [[nodiscard]] virtual double constraint_violation(double constraint_value, size_t constraint_index) const;
      template <typename Array>
      double constraint_violation(const Array& constraints, Norm residual_norm) const;
   };

   // compute ||c||
   template <typename Array>
   double Model::constraint_violation(const Array& constraints, Norm residual_norm) const {
      const Range constraints_range = Range(constraints.size());
      const VectorExpression constraint_violation{constraints_range, [&](size_t constraint_index) {
         return this->constraint_violation(constraints[constraint_index], constraint_index);
      }};
      return norm(residual_norm, constraint_violation);
   }
} // namespace

#endif // UNO_MODEL_H