// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"

namespace uno {
   // exact Hessian
   ExactHessian::ExactHessian(): HessianModel() {
   }

   void ExactHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const { }

   void ExactHessian::evaluate(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& primal_variables,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) {
      // evaluate Lagrangian Hessian
      hessian.set_dimension(problem.number_variables);
      problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, hessian);
      this->evaluation_count++;
   }
} // namespace