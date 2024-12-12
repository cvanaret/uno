// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <string>
#include "InequalityHandlingMethod.hpp"
#include "InequalityHandlingMethodFactory.hpp"
#include "inequality_constrained_methods/QPMethod.hpp"
#include "inequality_constrained_methods/LPMethod.hpp"
#include "interior_point_methods/PrimalDualInteriorPointMethod.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"

namespace uno {
   std::unique_ptr<InequalityHandlingMethod> InequalityHandlingMethodFactory::create(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) {
      const std::string inequality_handling_method = options.get_string("inequality_handling_method");
      // inequality-constrained methods
      if (inequality_handling_method == "QP") {
         return std::make_unique<QPMethod>(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
               number_hessian_nonzeros, options);
      }
      else if (inequality_handling_method == "LP") {
         return std::make_unique<LPMethod>(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
               options);
      }
      // interior-point method
      else if (inequality_handling_method == "primal_dual_interior_point") {
         return std::make_unique<PrimalDualInteriorPointMethod>(number_variables, number_constraints, number_jacobian_nonzeros,
               number_hessian_nonzeros, options);
      }
      throw std::invalid_argument("Subproblem strategy " + inequality_handling_method + " is not supported");
   }

   std::vector<std::string> InequalityHandlingMethodFactory::available_strategies() {
      std::vector<std::string> strategies{};
      if (not QPSolverFactory::available_solvers().empty()) {
         strategies.emplace_back("QP");
         strategies.emplace_back("LP");
      }
      else if (not LPSolverFactory::available_solvers().empty()) {
         strategies.emplace_back("LP");
      }
      if (not SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
         strategies.emplace_back("primal_dual_interior_point");
      }
      return strategies;
   }
} // namespace
