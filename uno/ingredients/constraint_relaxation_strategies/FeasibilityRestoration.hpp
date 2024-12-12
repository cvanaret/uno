// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "OptimalityProblem.hpp"
#include "l1RelaxedProblem.hpp"

namespace uno {
   enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

   class FeasibilityRestoration : public ConstraintRelaxationStrategy {
   public:
      FeasibilityRestoration(const Model& model, const Options& options);

      void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) override;

      [[nodiscard]] size_t maximum_number_variables() const override;
      [[nodiscard]] size_t maximum_number_constraints() const override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
            double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      // primal-dual residuals
      void compute_primal_dual_residuals(Iterate& iterate) override;
      void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;

   private:
      const OptimalityProblem optimality_problem;
      l1RelaxedProblem feasibility_problem;
      Phase current_phase{Phase::OPTIMALITY};
      const std::string& inequality_handling_method_name;
      const double linear_feasibility_tolerance;
      const bool switch_to_optimality_requires_linearized_feasibility;
      ProgressMeasures reference_optimality_progress{};
      Vector<double> reference_optimality_primals{};

      // delegating constructor
      FeasibilityRestoration(const Model& model, OptimalityProblem&& optimality_problem, l1RelaxedProblem&& feasibility_problem, const Options& options);

      [[nodiscard]] const OptimizationProblem& current_problem() const;
      void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            Direction& direction, WarmstartInformation& warmstart_information);
      void switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate, WarmstartInformation& warmstart_information);

      void evaluate_progress_measures(Iterate& iterate) const override;
      [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length);
      [[nodiscard]] bool can_switch_to_optimality_phase(const Iterate& current_iterate, const Iterate& trial_iterate, const Direction& direction,
            double step_length);
   };
} // namespace

#endif //UNO_FEASIBILITYRESTORATION_H
