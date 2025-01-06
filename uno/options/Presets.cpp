// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "Presets.hpp"
#include "Options.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"

namespace uno {
   Options Presets::get_preset_options(const std::optional<std::string>& optional_preset) {
      Options options(false);

      /** optional user preset **/
      if (optional_preset.has_value()) {
         Presets::set(options, *optional_preset);
      }
      else {
         /** default preset **/
         const auto QP_solvers = QPSolverFactory::available_solvers();
         const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
         const auto LP_solvers = LPSolverFactory::available_solvers();

         if (not QP_solvers.empty()) {
            Presets::set(options, "filtersqp");
         }
         else if (not linear_solvers.empty()) {
            Presets::set(options, "ipopt");
         }
         else if (not LP_solvers.empty()) {
            Presets::set(options, "filterslp");
         }
         // note: byrd is not very robust and is not considered as a default preset
      }
      return options;
   }

   void Presets::set(Options& options, const std::string& preset_name) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["inequality_handling_method"] = "primal_dual_interior_point";
         options["globalization_mechanism"] = "LS";
         options["globalization_strategy"] = "waechter_filter_method";
         options["filter_type"] = "standard";
         options["filter_beta"] = "0.99999";
         options["filter_gamma"] = "1e-8";
         options["switching_delta"] = "1";
         options["filter_ubd"] = "1e4";
         options["filter_fact"] = "1e4";
         options["filter_switching_infeasibility_exponent"] = "1.1";
         options["armijo_decrease_fraction"] = "1e-8";
         options["LS_backtracking_ratio"] = "0.5";
         options["LS_min_step_length"] = "5e-7";
         options["barrier_tau_min"] = "0.99";
         options["barrier_damping_factor"] = "1e-5";
         options["l1_constraint_violation_coefficient"] = "1000.";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "INF";
         options["scale_functions"] = "yes";
         options["sparse_format"] = "COO";
         options["tolerance"] = "1e-8";
         options["loose_tolerance"] = "1e-6";
         options["loose_tolerance_consecutive_iteration_threshold"] = "15";
         options["switch_to_optimality_requires_linearized_feasibility"] = "no";
         options["LS_scale_duals_with_step_length"] = "yes";
         options["protect_actual_reduction_against_roundoff"] = "yes";
      }
      else if (preset_name == "filtersqp") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["inequality_handling_method"] = "QP";
         options["globalization_mechanism"] = "TR";
         options["globalization_strategy"] = "fletcher_filter_method";
         options["filter_type"] = "standard";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L2";
         options["sparse_format"] = "CSC";
         options["TR_radius"] = "10";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["enforce_linear_constraints"] = "yes";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["TR_min_radius"] = "1e-8";
         options["switch_to_optimality_requires_linearized_feasibility"] = "yes";
         options["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "byrd") {
         options["constraint_relaxation_strategy"] = "l1_relaxation";
         options["inequality_handling_method"] = "QP";
         options["globalization_mechanism"] = "LS";
         options["globalization_strategy"] = "l1_merit";
         options["l1_relaxation_initial_parameter"] = "1";
         options["LS_backtracking_ratio"] = "0.5";
         options["armijo_decrease_fraction"] = "1e-8";
         options["l1_relaxation_epsilon1"] = "0.1";
         options["l1_relaxation_epsilon2"] = "0.1";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L1";
         options["sparse_format"] = "CSC";
         options["LS_scale_duals_with_step_length"] = "no";
         options["protect_actual_reduction_against_roundoff"] = "no";
      }
      else if (preset_name == "funnelsqp") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["inequality_handling_method"] = "QP";
         options["globalization_mechanism"] = "TR";
         options["globalization_strategy"] = "funnel_method";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L2";
         options["sparse_format"] = "CSC";
         options["TR_radius"] = "10";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["enforce_linear_constraints"] = "yes";
         options["tolerance"] = "1e-6";
         options["loose_tolerance"] = "1e-6";
         options["TR_min_radius"] = "1e-8";
         options["switch_to_optimality_requires_acceptance"] = "no";
         options["switch_to_optimality_requires_linearized_feasibility"] = "yes";

         options["funnel_beta"] = "0.9999";
         options["funnel_gamma"] = "0.001";
         options["switching_delta"] = "0.999";
         options["funnel_kappa"] = "0.5";
         options["funnel_ubd"] = "1.0";
         options["funnel_fact"] = "1.5";
         options["funnel_switching_infeasibility_exponent"] = "2";
         options["funnel_update_strategy"] = "2";
      }
      else if (preset_name == "filterslp") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["inequality_handling_method"] = "LP";
         options["globalization_mechanism"] = "TR";
         options["globalization_strategy"] = "fletcher_filter_method";
         options["filter_type"] = "standard";
         options["progress_norm"] = "L1";
         options["residual_norm"] = "L2";
         options["sparse_format"] = "CSC";
         options["TR_radius"] = "10";
         options["l1_constraint_violation_coefficient"] = "1.";
         options["enforce_linear_constraints"] = "yes";
         options["tolerance"] = "1e-5";
         options["loose_tolerance"] = "1e-4";
         options["TR_min_radius"] = "1e-8";
         options["switch_to_optimality_requires_linearized_feasibility"] = "yes";
         options["protect_actual_reduction_against_roundoff"] = "no";
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known");
      }
   }
} // namespace
