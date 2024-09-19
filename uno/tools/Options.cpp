// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <fstream>
#include <sstream>
#include "Options.hpp"
#include "Logger.hpp"

namespace uno {
   std::string& Options::operator[](const std::string& key) {
      return this->options[key];
   }

   const std::string& Options::at(const std::string& key) const {
      try {
         const std::string& value = this->options.at(key);
         this->is_used[key] = true;
         return value;
      }
      catch(const std::out_of_range&) {
         throw std::out_of_range("The option " + key + " was not found");
      }
   }

   const std::string& Options::get_string(const std::string& key) const {
      return this->at(key);
   }

   double Options::get_double(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stod(entry);
   }

   int Options::get_int(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stoi(entry);
   }

   size_t Options::get_unsigned_int(const std::string& key) const {
      const std::string& entry = this->at(key);
      return std::stoul(entry);
   }

   bool Options::get_bool(const std::string& key) const {
      const std::string& entry = this->at(key);
      return entry == "yes";
   }

   std::string Options::to_string(bool only_used) const {
      std::string representation = "Options:";
      for (const auto& [key, value]: this->options) {
         if (not only_used || this->is_used[key]) {
            representation += "\n- ";
            representation += key;
            representation += " = ";
            representation += value;
         }
      }
      return representation;
   }

   void Options::print(bool only_used) const {
      std::cout << this->to_string(only_used) << '\n';
   }

   Options Options::get_default_options(const std::string& file_name) {
      std::ifstream file;
      file.open(file_name);
      if (!file) {
         throw std::invalid_argument("The option file " + file_name + " was not found");
      }
      else {
         // register the default options
         Options options;
         std::string key, value;
         std::string line;
         while (std::getline(file, line)) {
            if (not line.empty() && line.find('#') != 0) {
               std::istringstream iss;
               iss.str(line);
               iss >> key >> value;
               options[key] = value;
            }
         }
         file.close();
         return options;
      }
   }

   void find_preset(const std::string& preset_name, Options& options) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         options["constraint_relaxation_strategy"] = "feasibility_restoration";
         options["subproblem"] = "primal_dual_interior_point";
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
         options["subproblem"] = "QP";
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
         options["subproblem"] = "QP";
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
         options["subproblem"] = "QP";
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
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known.");
      }
   }

   void Options::get_command_line_arguments(int argc, char* argv[]) {
      // build the (name, value) map
      int i = 1;
      while (i < argc - 1) {
         std::string argument = std::string(argv[i]);
         if (argument[0] == '-') {
            if (i < argc - 1) {
               // remove the '-'
               const std::string name = argument.substr(1);
               const std::string value = std::string(argv[i + 1]);
               if (name == "preset") {
                  find_preset(value, *this);
               }
               else {
                  this->operator[](name) = value;
               }
               i += 2;
            }
         }
         else {
            WARNING << "Argument " << argument << " was ignored\n";
            i++;
         }
      }
   }
} // namespace
