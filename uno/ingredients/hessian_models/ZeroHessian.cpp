// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"

namespace uno {
   void ZeroHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const { }

   void ZeroHessian::evaluate(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& /*primal_variables*/,
         const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(problem.number_variables);
      hessian.reset();
   }
}