// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QPMethod.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/QPSolver.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblems/LagrangeNewtonSubproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   QPMethod::QPMethod(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
         size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options) :
         InequalityConstrainedMethod(options.get_string("hessian_model"), number_variables, number_constraints, number_hessian_nonzeros,
               options.get_string("globalization_mechanism") != "TR" || options.get_bool("convexify_QP"), options),
         enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
         // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
         solver(QPSolverFactory::create(number_variables, number_constraints, number_objective_gradient_nonzeros, number_jacobian_nonzeros,
               // if the QP solver is used during preprocessing, we need to allocate the Hessian with at least number_variables elements
               std::max(this->enforce_linear_constraints_at_initial_iterate ? number_variables : 0, number_hessian_nonzeros),
               options)) {
   }

   QPMethod::~QPMethod() { }

   void QPMethod::generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) {
      if (this->enforce_linear_constraints_at_initial_iterate) {
         Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, *this->solver);
      }
   }

   void QPMethod::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
         Direction& direction, WarmstartInformation& warmstart_information) {
      LagrangeNewtonSubproblem subproblem(problem, current_iterate, current_multipliers.constraints, *this->hessian_model, this->trust_region_radius);
      this->solver->solve_QP(statistics, subproblem, this->initial_point, direction, warmstart_information);
      InequalityConstrainedMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
      // reset the initial point
      this->initial_point.fill(0.);
   }

   double QPMethod::hessian_quadratic_product(const Vector<double>& primal_direction) const {
      return this->solver->hessian_quadratic_product(primal_direction);
   }
} // namespace
