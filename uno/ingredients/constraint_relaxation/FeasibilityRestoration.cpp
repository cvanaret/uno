#include <cassert>
#include <cmath>
#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"

FeasibilityRestoration::FeasibilityRestoration(const Problem& problem, const Options& options) :
      ConstraintRelaxationStrategy(problem, options),
      // create the globalization strategies (one for each phase)
      phase_1_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)),
      phase_2_strategy(GlobalizationStrategyFactory::create(options.at("strategy"), options)) {
}

void FeasibilityRestoration::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("phase", Statistics::int_width, 4);

   // initialize the subproblem
   this->subproblem->initialize(statistics, problem, first_iterate);

   // compute the progress measures and the residuals of the initial point
   this->subproblem->compute_progress_measures(problem, first_iterate);
   this->subproblem->compute_residuals(problem, first_iterate, problem.objective_sign);

   // initialize the globalization strategies
   this->phase_1_strategy->initialize(statistics, first_iterate);
   this->phase_2_strategy->initialize(statistics, first_iterate);
}

void FeasibilityRestoration::create_current_subproblem(const Problem& problem, Iterate& current_iterate,
      double trust_region_radius) {
   this->subproblem->build_current_subproblem(problem, current_iterate, problem.objective_sign, trust_region_radius);
}

Direction FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, const Problem& problem, Iterate& current_iterate) {
   // solve the original subproblem
   Direction direction = this->subproblem->solve(statistics, problem, current_iterate);
   direction.objective_multiplier = problem.objective_sign;
   DEBUG << direction << "\n";
   assert((direction.status == OPTIMAL || direction.status == INFEASIBLE) && "The subproblem was not solved properly");

   // infeasible subproblem: form the feasibility problem
   if (direction.status == INFEASIBLE) {
      // try to minimize the constraint violation by solving the feasibility subproblem
      direction = this->solve_feasibility_problem(statistics, problem, current_iterate, direction.x, direction.constraint_partition);
      DEBUG << direction << "\n";
   }
   return direction;
}

Direction FeasibilityRestoration::solve_feasibility_problem(Statistics& statistics, const Problem& problem, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_primal_direction,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // form and solve the feasibility problem (with or without constraint partition)
   this->create_current_feasibility_problem(problem, current_iterate, optional_phase_2_primal_direction, optional_constraint_partition);

   DEBUG << "\nSolving the feasibility subproblem at the current iterate:\n" << current_iterate << "\n";
   Direction feasibility_direction = this->subproblem->solve(statistics, problem, current_iterate);
   feasibility_direction.objective_multiplier = 0.;
   DEBUG << feasibility_direction << "\n";
   assert(feasibility_direction.status == OPTIMAL && "The subproblem was not solved to optimality");

   if (optional_constraint_partition.has_value()) {
      // transfer the constraint partition of the phase-2 direction to the phase-1 direction
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      feasibility_direction.constraint_partition = constraint_partition;
   }
   else {
      // remove the temporary elastic variables
      this->remove_elastic_variables_from_subproblem();
   }
   return feasibility_direction;
}

bool FeasibilityRestoration::is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) {
   // check if subproblem definition changed
   if (this->subproblem->subproblem_definition_changed) {
      DEBUG << "The subproblem definition changed, the progress measures are recomputed\n";
      this->subproblem->subproblem_definition_changed = false;
      this->phase_2_strategy->reset();
      this->subproblem->compute_progress_measures(problem, current_iterate);
   }

   bool accept = false;
   if (ConstraintRelaxationStrategy::is_small_step(direction)) {
      this->subproblem->compute_progress_measures(problem, trial_iterate);
      accept = true;
   }
   else {
      // possibly switch between phase 1 (restoration) and phase 2 (optimality)
      GlobalizationStrategy& current_phase_strategy = this->switch_phase(problem, current_iterate, trial_iterate, direction);

      // evaluate the predicted reduction
      const double predicted_reduction = predicted_reduction_model.evaluate(step_length);

      // invoke the globalization strategy for acceptance
      accept = current_phase_strategy.check_acceptance(statistics, current_iterate.progress, trial_iterate.progress, direction.objective_multiplier,
            predicted_reduction);
   }

   if (accept) {
      statistics.add_statistic("phase", static_cast<int>(this->current_phase));
      // correct multipliers for infeasibility problem with constraint partition
      if (direction.objective_multiplier == 0. && direction.constraint_partition.has_value()) {
         const ConstraintPartition& constraint_partition = direction.constraint_partition.value();
         FeasibilityRestoration::set_restoration_multipliers(trial_iterate.multipliers.constraints, constraint_partition);
      }
      this->subproblem->compute_residuals(problem, trial_iterate, direction.objective_multiplier);
   }
   return accept;
}

void FeasibilityRestoration::add_proximal_term_to_subproblem(const Iterate& current_iterate) {
   // define a diagonal, inverse, quadratic proximal term
   this->subproblem->add_proximal_term_to_hessian([&](size_t i) {
      return std::pow(std::min(1., 1 / std::abs(current_iterate.x[i])), 2);
   });
}

void FeasibilityRestoration::create_current_feasibility_problem(const Problem& problem, Iterate& current_iterate,
      const std::optional<std::vector<double>>& optional_phase_2_primal_direction,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // if a constraint partition is given, form a partitioned l1 feasibility problem
   if (optional_constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      assert(!constraint_partition.infeasible.empty() && "The subproblem is infeasible but no constraint is infeasible");
      // set the multipliers of the violated constraints
      FeasibilityRestoration::set_restoration_multipliers(current_iterate.multipliers.constraints, constraint_partition);

      // compute the objective model with a zero objective multiplier
      this->subproblem->objective_gradient.clear();
      this->subproblem->build_objective_model(problem, current_iterate, 0.);
      if (this->use_proximal_term) {
         this->add_proximal_term_to_subproblem(current_iterate);
      }

      // assemble the linear objective (sum of the gradients of the violated constraints)
      this->subproblem->compute_feasibility_linear_objective(current_iterate, constraint_partition);

      // update the bounds of the constraints
      this->subproblem->generate_feasibility_bounds(problem, current_iterate.constraints, constraint_partition);
   }
   else {
      // no constraint partition given, form an l1 feasibility problem by adding elastic variables
      initialize_vector(current_iterate.multipliers.constraints, 0.);
      this->subproblem->build_objective_model(problem, current_iterate, 0.);
      if (this->use_proximal_term) {
         this->add_proximal_term_to_subproblem(current_iterate);
      }
      this->add_elastic_variables_to_subproblem(problem, current_iterate);
   }
   // start from the phase-2 solution
   if (optional_phase_2_primal_direction.has_value()) {
      const std::vector<double>& phase_2_primal_direction = optional_phase_2_primal_direction.value();
      this->subproblem->set_initial_point(phase_2_primal_direction);
   }
}

GlobalizationStrategy& FeasibilityRestoration::switch_phase(const Problem& problem, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction) {
   // possibly go from 1 (restoration) to phase 2 (optimality)
   if (this->current_phase == FEASIBILITY_RESTORATION && 0. < direction.objective_multiplier) {
      // TODO && this->filter_optimality->accept(trial_iterate.progress.feasibility, trial_iterate.progress.objective))
      this->current_phase = OPTIMALITY;
      DEBUG << "Switching from restoration to optimality phase\n";
      if (direction.constraint_partition.has_value()) {
         // remove elastics from the current iterate
         current_iterate.set_number_variables(this->number_subproblem_variables);
      }
      current_iterate.evaluate_constraints(problem);
      this->subproblem->compute_progress_measures(problem, current_iterate);
   }
   // possibly go from phase 2 (optimality) to 1 (restoration)
   else if (this->current_phase == OPTIMALITY && direction.objective_multiplier == 0.) {
      this->current_phase = FEASIBILITY_RESTORATION;
      DEBUG << "Switching from optimality to restoration phase\n";
      this->phase_2_strategy->notify(current_iterate);
      this->phase_1_strategy->reset();
      this->compute_infeasibility_measures(problem, current_iterate, direction.constraint_partition);
      this->phase_1_strategy->notify(current_iterate);
   }

   // evaluate the progress measures of the trial iterate
   if (this->current_phase == OPTIMALITY) {
      trial_iterate.set_number_variables(this->number_subproblem_variables);
      trial_iterate.evaluate_constraints(problem);
      this->subproblem->compute_progress_measures(problem, trial_iterate);
   }
   else { // restoration phase
      this->compute_infeasibility_measures(problem, trial_iterate, direction.constraint_partition);
      if (this->use_proximal_term) {
         this->add_proximal_term_to_progress_measures(current_iterate, trial_iterate);
      }
   }
   // return the globalization strategy of the current phase
   return (this->current_phase == OPTIMALITY) ? *this->phase_2_strategy : *this->phase_1_strategy;
}

void FeasibilityRestoration::set_restoration_multipliers(std::vector<double>& constraint_multipliers, const ConstraintPartition&
constraint_partition) {
   // the values {1, -1} are derived from the KKT conditions of the feasibility problem
   for (size_t j: constraint_partition.lower_bound_infeasible) {
      constraint_multipliers[j] = 1.;
   }
   for (size_t j: constraint_partition.upper_bound_infeasible) {
      constraint_multipliers[j] = -1.;
   }
   // otherwise, leave the multiplier as it is
}

void FeasibilityRestoration::compute_infeasibility_measures(const Problem& problem, Iterate& iterate,
      const std::optional<ConstraintPartition>& optional_constraint_partition) {
   // optimality measure: residual of linearly infeasible constraints
   if (optional_constraint_partition.has_value()) {
      const ConstraintPartition& constraint_partition = optional_constraint_partition.value();
      // feasibility measure: residual of all constraints
      iterate.evaluate_constraints(problem);
      const double feasibility_measure = problem.compute_constraint_violation(iterate.constraints, this->subproblem->residual_norm);
      const double objective_measure = problem.compute_constraint_violation(iterate.constraints, constraint_partition.infeasible,
            this->subproblem->residual_norm);
      iterate.progress = {feasibility_measure, objective_measure};
   }
   else {
      // if no constraint partition is available, simply compute the standard progress measures
      this->evaluate_relaxed_constraints(problem, iterate);
      this->subproblem->compute_progress_measures(problem, iterate);
      // add elastic variables to the optimality measure
      this->elastic_variables.negative.for_each_value([&](size_t i) {
         iterate.progress.objective += this->elastic_objective_coefficient*iterate.x[i];
      });
      this->elastic_variables.positive.for_each_value([&](size_t i) {
         iterate.progress.objective += this->elastic_objective_coefficient*iterate.x[i];
      });
   }
}

void FeasibilityRestoration::add_proximal_term_to_progress_measures(const Iterate& current_iterate, Iterate& trial_iterate) {
   const double coefficient = this->subproblem->get_proximal_coefficient();
   for (size_t i = 0; i < this->subproblem->number_variables; i++) {
      const double dr = std::min(1., 1/std::abs(current_iterate.x[i]));
      // measure weighted distance between trial iterate and current iterate
      const double proximal_term = coefficient * std::pow(dr*(trial_iterate.x[i] - current_iterate.x[i]), 2);
      trial_iterate.progress.objective += proximal_term;
   }
}