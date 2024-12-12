// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H
#define UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H

#include <memory>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // forward declaration
   class Options;

   // move the fixed variables to the set of general constraints (required in barrier methods)
   // for instance, the constraint x[2] == 1 is not interpreted as 1 <= x[2] <= 1, but as the
   // linear constraint 0 x[1] + 1 x[2] + ... + 0 x[n] == 1
   class FixedBoundsConstraintsModel: public Model {
   public:
      FixedBoundsConstraintsModel(std::unique_ptr<Model> original_model, const Options& options);

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override;
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override;
      void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
      void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const override;
      void compute_hessian_vector_product(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            Vector<double>& result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override;
      [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;

      void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const override;

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   private:
      const std::unique_ptr<Model> model{};
      Vector<size_t> fixed_variables;
      Vector<size_t> lower_bounded_variables;
      CollectionAdapter<Vector<size_t>&> lower_bounded_variables_collection;
      Vector<size_t> upper_bounded_variables;
      CollectionAdapter<Vector<size_t>&> upper_bounded_variables_collection;
      Concatenation<const Collection<size_t>&, ForwardRange> equality_constraints;
      Concatenation<const Collection<size_t>&, ForwardRange> linear_constraints;
   };
} // namespace

#endif // UNO_FIXEDBOUNDSCONSTRAINTSMODEL_H