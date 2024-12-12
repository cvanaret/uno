// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <functional>
#include "FeasibilityRestoration.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   FeasibilityRestoration::FeasibilityRestoration(const Model& model, const Options& options) :
         // call delegating constructor
         FeasibilityRestoration(model, OptimalityProblem(model),
               // create the (restoration phase) feasibility problem (objective multiplier = 0)
               l1RelaxedProblem(model, 0., options.get_double("l1_constraint_violation_coefficient"), 0., nullptr),
               options) {
   }

   // private delegating constructor
   FeasibilityRestoration::FeasibilityRestoration(const Model& model, OptimalityProblem&& optimality_problem, l1RelaxedProblem&& feasibility_problem,
            const Options& options) :
         ConstraintRelaxationStrategy(model,
               // allocate the largest size necessary to solve the optimality subproblem or the feasibility subproblem
               std::max(optimality_problem.number_variables, feasibility_problem.number_variables),
               std::max(optimality_problem.number_constraints, feasibility_problem.number_constraints),
               std::max(optimality_problem.number_objective_gradient_nonzeros(), feasibility_problem.number_objective_gradient_nonzeros()),
               std::max(optimality_problem.number_jacobian_nonzeros(), feasibility_problem.number_jacobian_nonzeros()),
               std::max(optimality_problem.number_hessian_nonzeros(), feasibility_problem.number_hessian_nonzeros()),
               options),
         optimality_problem(std::forward<OptimalityProblem>(optimality_problem)),
         feasibility_problem(std::forward<l1RelaxedProblem>(feasibility_problem)),
         inequality_handling_method_name(options.get_string("inequality_handling_method")),
         linear_feasibility_tolerance(options.get_double("tolerance")),
         switch_to_optimality_requires_linearized_feasibility(options.get_bool("switch_to_optimality_requires_linearized_feasibility")),
         reference_optimality_primals(optimality_problem.number_variables) {
      this->feasibility_problem.set_proximal_center(this->reference_optimality_primals.data());
   }

   void FeasibilityRestoration::initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) {
      // statistics
      this->inequality_handling_method->initialize_statistics(statistics, options);
      statistics.add_column("phase", Statistics::int_width, options.get_int("statistics_restoration_phase_column_order"));
      statistics.set("phase", "OPT");

      // initial iterate
      initial_iterate.feasibility_residuals.lagrangian_gradient.resize(this->feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.lower_bounds.resize(this->feasibility_problem.number_variables);
      initial_iterate.feasibility_multipliers.upper_bounds.resize(this->feasibility_problem.number_variables);
      this->inequality_handling_method->generate_initial_iterate(this->optimality_problem, initial_iterate);
      this->evaluate_progress_measures(initial_iterate);
      this->compute_primal_dual_residuals(initial_iterate);
      this->set_statistics(statistics, initial_iterate);
      this->globalization_strategy->initialize(statistics, initial_iterate, options);
   }

   void FeasibilityRestoration::compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) {
      direction.reset();
      // if we are in the optimality phase, solve the optimality problem
      if (this->current_phase == Phase::OPTIMALITY) {
         statistics.set("phase", "OPT");
         try {
            DEBUG << "Solving the optimality subproblem\n";
            this->solve_subproblem(statistics, this->optimality_problem, current_iterate, current_iterate.multipliers, direction, warmstart_information);
            if (direction.status == SubproblemStatus::INFEASIBLE) {
               // switch to the feasibility problem, starting from the current direction
               statistics.set("status", std::string("infeasible " + this->inequality_handling_method_name));
               DEBUG << "/!\\ The subproblem is infeasible\n";
               this->switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
               this->inequality_handling_method->set_initial_point(direction.primals);
            }
            else {
               warmstart_information.no_changes();
               return;
            }
         }
         catch (const UnstableRegularization&) {
            this->switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
         }
      }

      // solve the feasibility problem (minimize the constraint violation)
      DEBUG << "Solving the feasibility subproblem\n";
      statistics.set("phase", "FEAS");
      // note: failure of regularization should not happen here, since the feasibility Jacobian has full rank
      this->solve_subproblem(statistics, this->feasibility_problem, current_iterate, current_iterate.feasibility_multipliers, direction,
            warmstart_information);
      std::swap(direction.multipliers, direction.feasibility_multipliers);
   }

   bool FeasibilityRestoration::solving_feasibility_problem() const {
      return (this->current_phase == Phase::FEASIBILITY_RESTORATION);
   }

   // precondition: this->current_phase == Phase::OPTIMALITY
   void FeasibilityRestoration::switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) {
      DEBUG << "Switching from optimality to restoration phase\n";
      this->current_phase = Phase::FEASIBILITY_RESTORATION;
      this->globalization_strategy->notify_switch_to_feasibility(current_iterate.progress);
      this->inequality_handling_method->initialize_feasibility_problem(this->feasibility_problem, current_iterate);
      // save the current point (progress and primals) upon switching
      this->reference_optimality_progress = current_iterate.progress;
      this->reference_optimality_primals = current_iterate.primals;
      this->feasibility_problem.set_proximal_multiplier(this->inequality_handling_method->proximal_coefficient(current_iterate));

      current_iterate.set_number_variables(this->feasibility_problem.number_variables);
      this->inequality_handling_method->set_elastic_variable_values(this->feasibility_problem, current_iterate);
      DEBUG2 << "Current iterate:\n" << current_iterate << '\n';

      if (Logger::level == INFO) statistics.print_current_line();
      warmstart_information.whole_problem_changed();
   }

   void FeasibilityRestoration::solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, WarmstartInformation& warmstart_information) {
      direction.set_dimensions(problem.number_variables, problem.number_constraints);
      this->inequality_handling_method->solve(statistics, problem, current_iterate, current_multipliers, direction, warmstart_information);
      direction.norm = norm_inf(view(direction.primals, 0, this->model.number_variables));
      DEBUG3 << direction << '\n';
   }

   bool FeasibilityRestoration::can_switch_to_optimality_phase(const Iterate& current_iterate, const Iterate& trial_iterate, const Direction& direction,
         double step_length) {
      return this->globalization_strategy->is_infeasibility_sufficiently_reduced(this->reference_optimality_progress, trial_iterate.progress) &&
         (not this->switch_to_optimality_requires_linearized_feasibility ||
         this->model.constraint_violation(current_iterate.evaluations.constraints + step_length*(current_iterate.evaluations.constraint_jacobian *
         direction.primals), this->residual_norm) <= this->linear_feasibility_tolerance);
   }

   void FeasibilityRestoration::switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate, WarmstartInformation& warmstart_information) {
      DEBUG << "Switching from restoration to optimality phase\n";
      this->current_phase = Phase::OPTIMALITY;
      this->globalization_strategy->notify_switch_to_optimality(current_iterate.progress);
      current_iterate.set_number_variables(this->optimality_problem.number_variables);
      trial_iterate.set_number_variables(this->optimality_problem.number_variables);
      current_iterate.objective_multiplier = trial_iterate.objective_multiplier = 1.;

      this->inequality_handling_method->exit_feasibility_problem(this->optimality_problem, trial_iterate);
      // set a cold start in the subproblem solver
      warmstart_information.whole_problem_changed();
   }

   bool FeasibilityRestoration::is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      // TODO pick right multipliers
      this->inequality_handling_method->postprocess_iterate(this->current_problem(), trial_iterate);
      this->compute_progress_measures(current_iterate, trial_iterate);
      trial_iterate.objective_multiplier = this->current_problem().get_objective_multiplier();

      // possibly go from restoration phase to optimality phase
      if (this->current_phase == Phase::FEASIBILITY_RESTORATION && this->can_switch_to_optimality_phase(current_iterate, trial_iterate, direction, step_length)) {
         this->switch_to_optimality_phase(current_iterate, trial_iterate, warmstart_information);
      }
      else {
         warmstart_information.no_changes();
      }

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(this->model);
         accept_iterate = true;
         statistics.set("status", "0 primal step");
      }
      else {
         // invoke the globalization strategy for acceptance
         const ProgressMeasures predicted_reduction = this->compute_predicted_reduction_models(current_iterate, direction, step_length);
         accept_iterate = this->globalization_strategy->is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
               predicted_reduction, this->current_problem().get_objective_multiplier());
      }
      ConstraintRelaxationStrategy::set_progress_statistics(statistics, trial_iterate);
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals,
               this->current_phase == Phase::OPTIMALITY ? trial_iterate.multipliers : trial_iterate.feasibility_multipliers,
               this->current_problem().get_objective_multiplier());
      }
      return accept_iterate;
   }

   void FeasibilityRestoration::compute_primal_dual_residuals(Iterate& iterate) {
      ConstraintRelaxationStrategy::compute_primal_dual_residuals(this->optimality_problem, this->feasibility_problem, iterate);
   }

   const OptimizationProblem& FeasibilityRestoration::current_problem() const {
      if (this->current_phase == Phase::OPTIMALITY) {
         return this->optimality_problem;
      }
      else {
         return this->feasibility_problem;
      }
   }

   void FeasibilityRestoration::evaluate_progress_measures(Iterate& iterate) const {
      this->set_infeasibility_measure(iterate);
      this->set_objective_measure(iterate);
      this->inequality_handling_method->set_auxiliary_measure(this->model, iterate);
   }

   ProgressMeasures FeasibilityRestoration::compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length) {
      return {
         this->compute_predicted_infeasibility_reduction_model(current_iterate, direction.primals, step_length),
         this->compute_predicted_objective_reduction_model(current_iterate, direction.primals, step_length),
         this->inequality_handling_method->compute_predicted_auxiliary_reduction_model(this->model, current_iterate, direction.primals, step_length)
      };
   }

   size_t FeasibilityRestoration::maximum_number_variables() const {
      return std::max(this->optimality_problem.number_variables, this->feasibility_problem.number_variables);
   }

   size_t FeasibilityRestoration::maximum_number_constraints() const {
      return std::max(this->optimality_problem.number_constraints, this->feasibility_problem.number_constraints);
   }

   void FeasibilityRestoration::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const {
      const auto& residuals = (this->current_phase == Phase::OPTIMALITY) ? iterate.residuals : iterate.feasibility_residuals;
      statistics.set("stationarity", residuals.stationarity);
      statistics.set("complementarity", residuals.complementarity);
   }
} // namespace
