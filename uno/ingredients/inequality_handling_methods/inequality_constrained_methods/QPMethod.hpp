// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPMETHOD_H
#define UNO_QPMETHOD_H

#include <memory>
#include "InequalityConstrainedMethod.hpp"

namespace uno {
   // forward reference
   class QPSolver;

   class QPMethod : public InequalityConstrainedMethod {
   public:
      QPMethod(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, const Options& options);
      ~QPMethod();

      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
            Direction& direction, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& primal_direction) const override;

   protected:
      const bool enforce_linear_constraints_at_initial_iterate;
      // pointer to allow polymorphism
      const std::unique_ptr<QPSolver> solver;
   };
} // namespace

#endif // UNO_QPMETHOD_H