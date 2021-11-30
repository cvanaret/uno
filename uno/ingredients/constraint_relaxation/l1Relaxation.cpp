#include <cassert>
#include "l1Relaxation.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

l1Relaxation::l1Relaxation(Problem& problem, const Options& options) :
      ConstraintRelaxationStrategy(problem, options),
      globalization_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      penalty_parameter(stod(options.at("l1_relaxation_initial_parameter"))),
      parameters({stod(options.at("l1_relaxation_decrease_factor")),
                  stod(options.at("l1_relaxation_epsilon1")),
                  stod(options.at("l1_relaxation_epsilon2"))}),
      penalty_threshold(stod(options.at("l1_relaxation_penalty_threshold"))) {
}

void l1Relaxation::initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);

   // initialize the subproblem
   this->subproblem->initialize(statistics, problem, scaling, first_iterate);

   Subproblem::compute_optimality_conditions(problem, scaling, first_iterate, this->penalty_parameter);
   this->globalization_strategy->initialize(statistics, first_iterate);
}

void l1Relaxation::create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double trust_region_radius) {
   // scale the derivatives and introduce the elastic variables
   this->subproblem->create_current_subproblem(problem, scaling, current_iterate, this->penalty_parameter, trust_region_radius);
   this->add_elastic_variables_to_subproblem();

   // set the multipliers of the violated constraints
   l1Relaxation::set_multipliers(problem, current_iterate, this->subproblem->constraints_multipliers);
}

void l1Relaxation::set_multipliers(const Problem& problem, const Iterate& current_iterate, std::vector<double>& constraints_multipliers) {
   // the values are derived from the KKT conditions of the l1 problem
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (current_iterate.constraints[j] < problem.constraint_bounds[j].lb) { // lower bound infeasible
         constraints_multipliers[j] = 1.;
      }
      else if (problem.constraint_bounds[j].ub < current_iterate.constraints[j]) { // upper bound infeasible
         constraints_multipliers[j] = -1.;
      }
   }
   // otherwise, leave the multiplier as it is
}

Direction l1Relaxation::compute_feasible_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate) {
   DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";
   // use Byrd's steering rules to update the penalty parameter and compute descent directions
   Direction direction = this->solve_with_steering_rule(statistics, problem, scaling, current_iterate);

   // remove the temporary elastic variables from the direction
   this->remove_elastic_variables_from_direction(problem, direction);
   return direction;
}

Direction l1Relaxation::compute_second_order_correction(const Problem& problem, Iterate& trial_iterate) {
   Direction direction = ConstraintRelaxationStrategy::compute_second_order_correction(problem, trial_iterate);

   // remove the temporary elastic variables from the direction
   this->remove_elastic_variables_from_direction(problem, direction);
   return direction;
}

double l1Relaxation::compute_predicted_reduction(const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // compute the predicted reduction of the l1 relaxation as a postprocessing of the predicted reduction of the subproblem
   if (step_length == 1.) {
      return current_iterate.errors.constraints + predicted_reduction_model.evaluate(step_length);
   }
   else {
      // determine the linearized constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
      const auto residual_function = [&](size_t j) {
         const double component_j = current_iterate.constraints[j] + step_length * dot(direction.x, current_iterate.constraints_jacobian[j]);
         return problem.compute_constraint_violation(scaling, component_j, j);
      };
      const double linearized_constraint_violation = norm_1(residual_function, problem.number_constraints);
      return current_iterate.errors.constraints - linearized_constraint_violation + predicted_reduction_model.evaluate(step_length);
   }
}

Direction l1Relaxation::solve_feasibility_problem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      const Direction& /*phase_2_direction*/) {
   assert(0. < this->penalty_parameter && "l1Relaxation: the penalty parameter is already 0");

   const double objective_multiplier = 0.;
   Direction direction = this->resolve_subproblem(statistics, problem, scaling, current_iterate, objective_multiplier);
   // remove the temporary elastic variables
   this->remove_elastic_variables_from_direction(problem, direction);
   return direction;
}

bool l1Relaxation::is_acceptable(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      this->globalization_strategy->reset();
      this->subproblem->subproblem_definition_changed = false;
      this->subproblem->compute_progress_measures(problem, scaling, current_iterate);
   }

   bool accept = false;
   if (direction.norm == 0.) {
      accept = true;
   }
   else {
      this->subproblem->compute_progress_measures(problem, scaling, trial_iterate);

      // compute the predicted reduction (both the subproblem and the l1 relaxation strategy contribute)
      const double predicted_reduction = this->compute_predicted_reduction(problem, scaling, current_iterate, direction, predicted_reduction_model,
            step_length);
      // invoke the globalization strategy for acceptance
      accept = this->globalization_strategy->check_acceptance(statistics, current_iterate.progress, trial_iterate.progress,
            this->penalty_parameter, predicted_reduction);
   }
   if (accept) {
      statistics.add_statistic("penalty param.", this->penalty_parameter);
      Subproblem::compute_optimality_conditions(problem, scaling, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

Direction l1Relaxation::solve_with_steering_rule(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate) {
   // stage a: compute the step within trust region
   Direction direction = this->solve_subproblem(statistics, problem, current_iterate);

   // penalty update: if penalty parameter is already 0, no need to decrease it
   if (0. < this->penalty_parameter) {
      // check infeasibility
      double linearized_residual = this->compute_linearized_constraint_residual(direction.x);
      DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

      // if problem had to be relaxed
      if (linearized_residual != 0.) {
         const double current_penalty_parameter = this->penalty_parameter;

         // stage c: compute the lowest possible constraint violation (penalty = 0)
         DEBUG << "Compute ideal solution (param = 0):\n";
         Direction direction_lowest_violation = this->resolve_subproblem(statistics, problem, scaling, current_iterate, 0.);
         const double residual_lowest_violation = this->compute_linearized_constraint_residual(direction_lowest_violation.x);
         DEBUG << "Ideal linearized residual mk(dk): " << residual_lowest_violation << "\n\n";

         if (!(0. < current_iterate.errors.constraints && residual_lowest_violation == current_iterate.errors.constraints)) {
            // compute the ideal error (with a zero penalty parameter)
            const double error_lowest_violation = l1Relaxation::compute_error(problem, scaling, current_iterate,
                  direction_lowest_violation.multipliers, 0.);
            DEBUG << "Ideal error: " << error_lowest_violation << "\n";
            if (error_lowest_violation == 0.) {
               // stage f: update the penalty parameter
               this->penalty_parameter = 0.;
               direction = direction_lowest_violation;
            }
            else {
               // stage f: update the penalty parameter
               const double updated_penalty_parameter = this->penalty_parameter;
               const double term = error_lowest_violation / std::max(1., current_iterate.errors.constraints);
               this->penalty_parameter = std::min(this->penalty_parameter, term * term);
               if (this->penalty_parameter < updated_penalty_parameter) {
                  if (this->penalty_parameter == 0.) {
                     direction = direction_lowest_violation;
                  }
                  else {
                     direction = this->resolve_subproblem(statistics, problem, scaling, current_iterate, this->penalty_parameter);
                  }
               }

               // decrease penalty parameter to satisfy 2 conditions
               bool condition1 = false, condition2 = false;
               while (!condition2) {
                  if (!condition1) {
                     // stage d: reach a fraction of the ideal decrease
                     if ((residual_lowest_violation == 0. && linearized_residual == 0) || (residual_lowest_violation != 0. &&
                     current_iterate.errors.constraints - linearized_residual >= this->parameters.epsilon1 *
                     (current_iterate.errors.constraints - residual_lowest_violation))) {
                        condition1 = true;
                        DEBUG << "Condition 1 is true\n";
                     }
                  }
                  // stage e: further decrease penalty parameter if necessary
                  if (condition1 && current_iterate.errors.constraints - direction.objective >=
                                    this->parameters.epsilon2 * (current_iterate.errors.constraints - direction_lowest_violation.objective)) {
                     condition2 = true;
                     DEBUG << "Condition 2 is true\n";
                  }
                  if (!condition2) {
                     this->penalty_parameter /= this->parameters.decrease_factor;
                     if (this->penalty_parameter < this->penalty_threshold) {
                        this->penalty_parameter = 0.;
                        condition2 = true;
                     }
                     else {
                        DEBUG << "\nAttempting to solve with penalty parameter " << this->penalty_parameter << "\n";
                        direction = this->resolve_subproblem(statistics, problem, scaling, current_iterate, this->penalty_parameter);
   
                        linearized_residual = this->compute_linearized_constraint_residual(direction.x);
                        DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                     }
                  }
               }
            } // end else
         }

         if (this->penalty_parameter < current_penalty_parameter) {
            DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
            this->globalization_strategy->reset();
         }
      }
   }
   return direction;
}

Direction l1Relaxation::solve_subproblem(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   if (direction.constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      assert(constraint_partition.infeasible.empty() && "solve_subproblem: infeasible constraints found, although direction is feasible");
   }
   direction.objective_multiplier = this->penalty_parameter;
   DEBUG << "\n" << direction;

   // remove the temporary elastic variables
   this->remove_elastic_variables_from_subproblem();
   return direction;
}

Direction l1Relaxation::resolve_subproblem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
      double objective_multiplier) {
   this->subproblem->build_objective_model(problem, scaling, current_iterate, objective_multiplier);
   this->add_elastic_variables_to_subproblem();

   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   if (direction.constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      assert(constraint_partition.infeasible.empty() && "resolve_subproblem: infeasible constraints found, although direction is feasible");
   }
   direction.objective_multiplier = objective_multiplier;
   DEBUG << "\n" << direction;

   // remove the temporary elastic variables
   this->remove_elastic_variables_from_subproblem();
   return direction;
}

double l1Relaxation::compute_linearized_constraint_residual(std::vector<double>& direction) const {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   auto add_variable_contribution = [&](size_t i) {
      residual += direction[i];
   };
   this->elastic_variables.positive.for_each_value(add_variable_contribution);
   this->elastic_variables.negative.for_each_value(add_variable_contribution);
   return residual;
}

// measure that combines KKT error and complementarity error
double l1Relaxation::compute_error(const Problem& problem, const Scaling& scaling, Iterate& iterate, Multipliers& multipliers_displacements, double
current_penalty_parameter) {
   Multipliers multipliers(multipliers_displacements.lower_bounds.size(), multipliers_displacements.constraints.size());
   multipliers.lower_bounds = multipliers_displacements.lower_bounds;
   multipliers.upper_bounds = multipliers_displacements.upper_bounds;
   multipliers.constraints = multipliers_displacements.constraints;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      multipliers.constraints[j] += iterate.multipliers.constraints[j];
   }

   // complementarity error
   double error = Subproblem::compute_complementarity_error(problem, scaling, iterate, multipliers);
   // KKT error
   iterate.evaluate_lagrangian_gradient(problem, scaling, current_penalty_parameter, multipliers);
   error += norm_1(iterate.lagrangian_gradient);
   return error;
}