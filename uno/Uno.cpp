// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanism.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "optimization/OptimizationStatus.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
         globalization_mechanism(globalization_mechanism),
         max_iterations(options.get_unsigned_int("max_iterations")),
         time_limit(options.get_double("time_limit")),
         print_solution(options.get_bool("print_solution")),
         strategy_combination(Uno::get_strategy_combination(options)) { }
   
   Level Logger::level = INFO;

   // solve without user callbacks
   Result Uno::solve(const Model& model, Iterate& current_iterate, const Options& options) {
      // pass user callbacks that do nothing
      NoUserCallbacks user_callbacks{};
      return this->solve(model, current_iterate, options, user_callbacks);
   }

   // solve with user callbacks
   Result Uno::solve(const Model& model, Iterate& current_iterate, const Options& options, UserCallbacks& user_callbacks) {
      Timer timer{};
      Statistics statistics = Uno::create_statistics(model, options);
      WarmstartInformation warmstart_information{};
      warmstart_information.whole_problem_changed();

      size_t major_iterations = 0;
      OptimizationStatus optimization_status = OptimizationStatus::SUCCESS;
      try {
         // use the initial primal-dual point to initialize the strategies and generate the initial iterate
         this->initialize(statistics, current_iterate, options);
         // allocate the trial iterate once and for all here
         Iterate trial_iterate(current_iterate);

         try {
            bool termination = false;
            // check for termination
            while (not termination) {
               major_iterations++;
               statistics.start_new_line();
               statistics.set("iter", major_iterations);
               DEBUG << "### Outer iteration " << major_iterations << '\n';

               // compute an acceptable iterate by solving a subproblem at the current point
               warmstart_information.iterate_changed();
               this->globalization_mechanism.compute_next_iterate(statistics, model, current_iterate, trial_iterate, warmstart_information, user_callbacks);
               termination = this->termination_criteria(trial_iterate.status, major_iterations, timer.get_duration(), optimization_status);
               user_callbacks.notify_new_primals(trial_iterate.primals);
               user_callbacks.notify_new_multipliers(trial_iterate.multipliers);

               // the trial iterate becomes the current iterate for the next iteration
               std::swap(current_iterate, trial_iterate);
            }
         }
         catch (std::exception& exception) {
            statistics.start_new_line();
            statistics.set("status", exception.what());
            if (Logger::level == INFO) statistics.print_current_line();
            DEBUG << exception.what() << '\n';
            optimization_status = OptimizationStatus::ALGORITHMIC_ERROR;
         }
         if (Logger::level == INFO) statistics.print_footer();

         Uno::postprocess_iterate(model, current_iterate, current_iterate.status);
      }
      catch (const std::exception& e) {
         DISCRETE  << "An error occurred at the initial iterate: " << e.what()  << '\n';
         optimization_status = OptimizationStatus::EVALUATION_ERROR;
      }
      Result result = this->create_result(model, optimization_status, current_iterate, major_iterations, timer);
      this->print_optimization_summary(result);
      return result;
   }

   void Uno::initialize(Statistics& statistics, Iterate& current_iterate, const Options& options) {
      statistics.start_new_line();
      statistics.set("iter", 0);
      statistics.set("status", "initial point");
      this->globalization_mechanism.initialize(statistics, current_iterate, options);
      options.print_used();
      if (Logger::level == INFO) statistics.print_current_line();
      current_iterate.status = IterateStatus::NOT_OPTIMAL;
   }

   Statistics Uno::create_statistics(const Model& model, const Options& options) {
      Statistics statistics(options);
      statistics.add_column("iter", Statistics::int_width, options.get_int("statistics_major_column_order"));
      statistics.add_column("step norm", Statistics::double_width - 5, options.get_int("statistics_step_norm_column_order"));
      statistics.add_column("objective", Statistics::double_width - 5, options.get_int("statistics_objective_column_order"));
      if (model.is_constrained()) {
         statistics.add_column("primal feas", Statistics::double_width - 4, options.get_int("statistics_primal_feasibility_column_order"));
      }
      statistics.add_column("stationarity", Statistics::double_width - 3, options.get_int("statistics_stationarity_column_order"));
      statistics.add_column("complementarity", Statistics::double_width, options.get_int("statistics_complementarity_column_order"));
      statistics.add_column("status", Statistics::string_width - 9, options.get_int("statistics_status_column_order"));
      return statistics;
   }

   bool Uno::termination_criteria(IterateStatus current_status, size_t iteration, double current_time, OptimizationStatus& optimization_status) const {
      if (current_status != IterateStatus::NOT_OPTIMAL) {
         return true;
      }
      else if (this->max_iterations <= iteration) {
         optimization_status = OptimizationStatus::ITERATION_LIMIT;
         return true;
      }
      else if (this->time_limit <= current_time) {
         optimization_status = OptimizationStatus::TIME_LIMIT;
         return true;
      }
      return false;
   }

   void Uno::postprocess_iterate(const Model& model, Iterate& iterate, IterateStatus termination_status) {
      // in case the objective was not yet evaluated, evaluate it
      iterate.evaluate_objective(model);
      model.postprocess_solution(iterate, termination_status);
      DEBUG2 << "Final iterate:\n" << iterate;
   }

   Result Uno::create_result(const Model& model, OptimizationStatus optimization_status, Iterate& current_iterate, size_t major_iterations,
         const Timer& timer) {
      const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
      const size_t number_hessian_evaluations = this->globalization_mechanism.get_hessian_evaluation_count();
      return {optimization_status, std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations,
            timer.get_duration(), Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_objective_gradient,
            Iterate::number_eval_jacobian, number_hessian_evaluations, number_subproblems_solved};
   }

   std::string Uno::current_version() {
      return "1.3.0";
   }

   void Uno::print_available_strategies() {
      std::cout << "Available strategies:\n";
      std::cout << "- Constraint relaxation strategies: " << join(ConstraintRelaxationStrategyFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Globalization mechanisms: " << join(GlobalizationMechanismFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Globalization strategies: " << join(GlobalizationStrategyFactory::available_strategies(), ", ") << '\n';
      std::cout << "- Subproblems: " << join(InequalityHandlingMethodFactory::available_strategies(), ", ") << '\n';
      std::cout << "- QP solvers: " << join(QPSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- LP solvers: " << join(LPSolverFactory::available_solvers(), ", ") << '\n';
      std::cout << "- Linear solvers: " << join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", ") << '\n';
   }

   std::string Uno::get_strategy_combination(const Options& options) {
      return options.get_string("globalization_mechanism") + " " + options.get_string("constraint_relaxation_strategy") + " " +
                                options.get_string("globalization_strategy") + " " + options.get_string("inequality_handling_method");

   }

   void Uno::print_optimization_summary(const Result& result) {
      DISCRETE << "\nUno " << Uno::current_version() << " (" << this->strategy_combination << ")\n";
      DISCRETE << Timer::get_current_date();
      DISCRETE << "────────────────────────────────────────\n";
      result.print(this->print_solution);
   }
} // namespace
