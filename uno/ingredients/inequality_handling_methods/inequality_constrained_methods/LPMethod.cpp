// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPMethod.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "ingredients/subproblem_solvers/LPSolver.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"

namespace uno {
   LPMethod::LPMethod(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, const Options& options) :
         InequalityConstrainedMethod("zero", number_variables, number_constraints, 0, false, options),
         solver(LPSolverFactory::create(number_variables, number_constraints,
               number_objective_gradient_nonzeros, number_jacobian_nonzeros, options)) {
   }

   LPMethod::~LPMethod() { }

   void LPMethod::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
   }

   void LPMethod::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, WarmstartInformation& warmstart_information) {
      LagrangeNewtonSubproblem subproblem(problem, current_iterate, current_multipliers.constraints, *this->hessian_model, this->trust_region_radius);
      this->solver->solve_LP(statistics, subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   double LPMethod::hessian_quadratic_product(const Vector<double>& /*primal_direction*/) const {
      return 0.;
   }
} // namespace