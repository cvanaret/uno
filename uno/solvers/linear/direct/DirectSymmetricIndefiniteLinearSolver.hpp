// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H
#define UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H

#include <stdexcept>
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   template <typename NumericalType>
   using LinearOperator = void(const Vector<NumericalType>&, Vector<NumericalType>&);

   template <typename IndexType, typename NumericalType>
   class DirectSymmetricIndefiniteLinearSolver: public SymmetricIndefiniteLinearSolver<IndexType, NumericalType, LinearOperator<NumericalType>> {
   public:
      explicit DirectSymmetricIndefiniteLinearSolver(size_t dimension) :
         SymmetricIndefiniteLinearSolver<IndexType, NumericalType, LinearOperator<NumericalType>>(dimension) { };
      ~DirectSymmetricIndefiniteLinearSolver() = default;

      virtual void factorize(const SymmetricMatrix<IndexType, NumericalType>& matrix) = 0;
      virtual void do_symbolic_factorization(const SymmetricMatrix<IndexType, NumericalType>& matrix) = 0;
      virtual void do_numerical_factorization(const SymmetricMatrix<IndexType, NumericalType>& matrix) = 0;

      virtual void solve_indefinite_system(const SymmetricMatrix<IndexType, NumericalType>& matrix, const Vector<NumericalType>& rhs,
            Vector<NumericalType>& result) = 0;
      void solve_indefinite_system(const LinearOperator<NumericalType>& /*linear_operator*/, const Vector<NumericalType>& /*rhs*/,
            Vector<NumericalType>& /*result*/) override {
         throw std::runtime_error("DirectSymmetricIndefiniteLinearSolver: solve_indefinite_system with linear operator is not implemented yet.");
      }

      [[nodiscard]] virtual std::tuple<size_t, size_t, size_t> get_inertia() const = 0;
      [[nodiscard]] virtual size_t number_negative_eigenvalues() const = 0;
      // [[nodiscard]] virtual bool matrix_is_positive_definite() const = 0;
      [[nodiscard]] virtual bool matrix_is_singular() const = 0;
      [[nodiscard]] virtual size_t rank() const = 0;
   };
} // namespace

#endif // UNO_DIRECTSYMMETRICINDEFINITELINEARSOLVER_H