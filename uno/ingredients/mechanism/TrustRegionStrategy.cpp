#include <cmath>
#include <cassert>
#include "TrustRegionStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Logger.hpp"

TrustRegionStrategy::TrustRegionStrategy(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) :
      GlobalizationMechanism(constraint_relaxation_strategy),
      radius(stod(options.at("TR_radius"))),
      increase_factor(stod(options.at("TR_increase_factor"))),
      decrease_factor(stod(options.at("TR_decrease_factor"))),
      activity_tolerance(stod(options.at("TR_activity_tolerance"))),
      min_radius(stod(options.at("TR_min_radius"))) {
}

void TrustRegionStrategy::initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) {
   statistics.add_column("TR radius", Statistics::double_width, 30);
   // generate the initial point
   this->relaxation_strategy.initialize(statistics, problem, scaling, first_iterate);
}

std::tuple<Iterate, double> TrustRegionStrategy::compute_acceptable_iterate(Statistics& statistics, const Problem& problem, const Scaling& scaling,
      Iterate& current_iterate) {
   this->number_iterations = 0;

   while (!this->termination()) {
      try {
         assert (0 < this->radius);
         this->number_iterations++;
         this->print_iteration();

         // generate the subproblem
         this->relaxation_strategy.create_current_subproblem(problem, scaling, current_iterate, this->radius);

         // compute the direction within the trust region
         Direction direction = this->relaxation_strategy.compute_feasible_direction(statistics, problem, scaling, current_iterate);
         TrustRegionStrategy::check_unboundedness(direction);
         // set bound multipliers of active trust region to 0
         TrustRegionStrategy::rectify_active_set(direction, this->radius);

         // assemble the trial iterate by taking a full step
         const double full_step_length = 1.;
         Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, full_step_length);

         // check whether the trial step is accepted
         PredictedReductionModel predicted_reduction_model = this->relaxation_strategy.generate_predicted_reduction_model(problem, direction);
         if (this->relaxation_strategy.is_acceptable(statistics, problem, scaling, current_iterate, trial_iterate, direction,
               predicted_reduction_model, full_step_length)) {
            this->add_statistics(statistics, direction);

            // increase the radius if trust region is active
            if (direction.norm >= this->radius - this->activity_tolerance) {
               this->radius *= this->increase_factor;
            }

            // let the subproblem know the accepted iterate
            this->relaxation_strategy.register_accepted_iterate(trial_iterate);
            return std::make_tuple(std::move(trial_iterate), direction.norm);
         }
         else {
            // if the step is rejected, decrease the radius
            this->radius = std::min(this->radius, direction.norm) / this->decrease_factor;
         }
      }
      catch (const NumericalError& e) {
         GlobalizationMechanism::print_warning(e.what());
         // if an evaluation error occurs, decrease the radius
         this->radius /= this->decrease_factor;
      }
   }
   // radius gets too small
   if (this->radius < this->min_radius) { // 1e-16: something like machine precision
      throw std::runtime_error("Trust-region radius became too small");
   }
   else {
      throw std::runtime_error("Trust-region failed with an unexpected error");
   }
}

void TrustRegionStrategy::check_unboundedness(const Direction& direction) {
   assert(direction.status != UNBOUNDED_PROBLEM && "Trust-region subproblem is unbounded, this should not happen");
}

void TrustRegionStrategy::rectify_active_set(Direction& direction, double radius) {
   assert(0 < radius);
   // update active set and set multipliers for bound constraints active at trust region to 0
   for (auto it = direction.active_set.bounds.at_lower_bound.begin(); it != direction.active_set.bounds.at_lower_bound.end();) {
      size_t i = *it;
      if (direction.x[i] == -radius) {
         it = direction.active_set.bounds.at_lower_bound.erase(it);
         direction.multipliers.lower_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
   for (auto it = direction.active_set.bounds.at_upper_bound.begin(); it != direction.active_set.bounds.at_upper_bound.end();) {
      size_t i = *it;
      if (direction.x[i] == radius) {
         it = direction.active_set.bounds.at_upper_bound.erase(it);
         direction.multipliers.upper_bounds[i] = 0.;
      }
      else {
         ++it;
      }
   }
}

void TrustRegionStrategy::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("TR radius", this->radius);
   statistics.add_statistic("step norm", direction.norm);
}

bool TrustRegionStrategy::termination() {
   return this->radius < this->min_radius;
}

void TrustRegionStrategy::print_iteration() {
   DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
}