#include "QPSubproblem.hpp"
#include "solvers/QP/QPSolverFactory.hpp"

QPSubproblem::QPSubproblem(const NonlinearReformulation& problem, const Options& options) :
      ActiveSetSubproblem(problem, NO_SOC),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(options.at("QP_solver"), problem.number_variables, problem.number_constraints,
            problem.get_maximum_number_hessian_nonzeros()
            + problem.number_variables, /* regularization */
            true)),
      proximal_coefficient(stod(options.at("proximal_coefficient"))),
      // if no trust region is used, the problem should be convexified to guarantee boundedness + descent direction
      hessian_model(HessianModelFactory::create(options.at("hessian_model"), problem.number_variables,
            problem.get_maximum_number_hessian_nonzeros() + problem.number_variables, options.at("mechanism") != "TR", options)) {
}

void QPSubproblem::evaluate_problem(const NonlinearReformulation& problem, Iterate& current_iterate) {
   // Hessian
   this->hessian_model->evaluate(problem, current_iterate.primals, current_iterate.multipliers.constraints);

   // objective gradient
   problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);

   // constraints
   problem.evaluate_constraints(current_iterate, this->constraints);

   // constraint Jacobian
   problem.evaluate_constraint_jacobian(current_iterate, this->constraint_jacobian);
}

Direction QPSubproblem::solve(Statistics& /*statistics*/, const NonlinearReformulation& problem, Iterate& current_iterate) {
   // evaluate the functions at the current iterate
   this->evaluate_problem(problem, current_iterate);

   // bounds of the variable displacements
   this->set_variable_displacement_bounds(problem, current_iterate);

   // bounds of the linearized constraints
   this->set_linearized_constraint_bounds(problem, this->constraints);

   // compute QP direction
   Direction direction = this->solver->solve_QP(problem.number_variables, problem.number_constraints, this->variable_displacement_bounds,
         this->linearized_constraint_bounds, this->objective_gradient, this->constraint_jacobian, *this->hessian_model->hessian, this->initial_point);
   ActiveSetSubproblem::compute_dual_displacements(problem, current_iterate, direction);
   this->number_subproblems_solved++;
   return direction;
}

PredictedReductionModel QPSubproblem::generate_predicted_reduction_model(const NonlinearReformulation& problem, const Direction& direction) const {
   return PredictedReductionModel(-direction.objective, [&]() { // capture "this" and "direction" by reference
      // precompute expensive quantities
      const double linear_term = dot(direction.primals, this->objective_gradient);
      const double quadratic_term = this->hessian_model->hessian->quadratic_product(direction.primals, direction.primals, problem.number_variables) / 2.;
      // return a function of the step length that cheaply assembles the predicted reduction
      return [=](double step_length) { // capture the expensive quantities by value
         return -step_length * (linear_term + step_length * quadratic_term);
      };
   });
}

size_t QPSubproblem::get_hessian_evaluation_count() const {
   return this->hessian_model->evaluation_count;
}

double QPSubproblem::get_proximal_coefficient() const {
   return this->proximal_coefficient;
}