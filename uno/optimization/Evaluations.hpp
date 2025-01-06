// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONS_H
#define UNO_EVALUATIONS_H

#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   struct Evaluations {
      double objective{INF<double>}; /*!< Objective value */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */

      Evaluations(size_t number_variables, size_t number_constraints):
            constraints(number_constraints),
            objective_gradient(number_variables),
            constraint_jacobian(number_constraints, number_variables) {
      }
   };
} // namespace

#endif // UNO_EVALUATIONS_H
