#include <cassert>
#include "SQP.hpp"
#include "QPSolverFactory.hpp"

SQP::SQP(size_t number_variables, size_t number_constraints, size_t hessian_maximum_number_nonzeros, const std::string& QP_solver_name,
      const std::string& hessian_evaluation_method, bool use_trust_region) :
      Subproblem(number_variables, number_constraints),
      // maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
      solver(QPSolverFactory::create(QP_solver_name, number_variables, number_constraints,
            hessian_maximum_number_nonzeros + number_variables, true)),
      /* if no trust region is used, the problem should be convexified by controlling the inertia of the Hessian */
      hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, number_variables, hessian_maximum_number_nonzeros, !use_trust_region)),
      initial_point(number_variables) {
}

void SQP::generate(const Problem& problem, Iterate& current_iterate, double objective_multiplier, double trust_region_radius) {
   copy_from(this->constraints_multipliers, current_iterate.multipliers.constraints);
   /* compute first- and second-order information */
   problem.evaluate_constraints(current_iterate.x, current_iterate.constraints);
   this->constraints_jacobian = problem.constraints_jacobian(current_iterate.x);

   this->objective_gradient = problem.objective_gradient(current_iterate.x);
   this->update_objective_multiplier(problem, current_iterate, objective_multiplier);

   /* bounds of the variables */
   this->set_variables_bounds(problem, current_iterate, trust_region_radius);

   /* bounds of the linearized constraints */
   this->set_constraints_bounds(problem, current_iterate.constraints);

   /* set the initial point */
   clear(this->initial_point);
}

void SQP::update_objective_multiplier(const Problem& problem, const Iterate& current_iterate, double objective_multiplier) {
   // evaluate the Hessian
   this->hessian_evaluation->compute(problem, current_iterate.x, objective_multiplier, this->constraints_multipliers);

   // scale objective gradient
   if (objective_multiplier == 0.) {
      clear(this->objective_gradient);
   }
   else if (objective_multiplier < 1.) {
      this->objective_gradient = problem.objective_gradient(current_iterate.x);
      scale(this->objective_gradient, objective_multiplier);
   }
   clear(this->initial_point);
}

void SQP::set_initial_point(const std::vector<double>& point) {
   copy_from(this->initial_point, point);
}

Direction SQP::compute_direction(Statistics& /*statistics*/, const Problem& /*problem*/, Iterate& /*current_iterate*/) {
   /*
   std::cout << "Hessian:\n" << this->hessian_evaluation->hessian << "\n";
   std::cout << "Obj gradient: "; print_vector(std::cout, this->objective_gradient);
   std::cout << "Variables bounds:\n";
   for (const Range range: variables_bounds) {
      std::cout << range.lb << ", " << range.ub << "\n";
   }
   std::cout << "Constraints bounds:\n";
   for (const Range range: constraints_bounds) {
      std::cout << range.lb << ", " << range.ub << "\n";
   }
   for (size_t j = 0; j < this->constraints_jacobian.size(); j++) {
      std::cout << "Constraint " << j << " gradient: "; print_vector(std::cout, this->constraints_jacobian[j]);
   }
   std::cout << "Initial point: "; print_vector(std::cout, this->initial_point);
    */

   /* compute QP direction */
   Direction direction = this->solver->solve_QP(this->variables_bounds, this->constraints_bounds, this->objective_gradient,
         this->constraints_jacobian, this->hessian_evaluation->hessian, this->initial_point);
   this->number_subproblems_solved++;
   DEBUG << direction;

   // attach the predicted reduction function
   direction.predicted_reduction = [&](double step_length) {
      return this->compute_predicted_reduction(direction, step_length);
   };
   return direction;
}

double SQP::compute_predicted_reduction(const Direction& direction, double step_length) const {
   // the predicted reduction is quadratic in the step length
   if (step_length == 1.) {
      return -direction.objective;
   }
   else {
      double linear_term = dot(direction.x, this->objective_gradient);
      double quadratic_term = this->hessian_evaluation->hessian.quadratic_product(direction.x, direction.x) / 2.;
      return -step_length * (linear_term + step_length * quadratic_term);
   }
}

int SQP::get_hessian_evaluation_count() const {
   return this->hessian_evaluation->evaluation_count;
}