// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "AMPLModel.hpp"
#include "Uno.hpp"
#include "model/ModelFactory.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"

/*
size_t memory_allocation_amount = 0;

void* operator new(size_t size) {
   memory_allocation_amount += size;
   std::cout << "Memory: " << size << '\n';
   return malloc(size);
}
*/

namespace uno {
   void run_uno_ampl(const std::string& model_name, const Options& options) {
      // AMPL model
      std::unique_ptr<Model> ampl_model = std::make_unique<AMPLModel>(model_name);

      // initialize initial primal and dual points
      Iterate initial_iterate(ampl_model->number_variables, ampl_model->number_constraints);
      ampl_model->initial_primal_point(initial_iterate.primals);
      ampl_model->project_onto_variable_bounds(initial_iterate.primals);
      ampl_model->initial_dual_point(initial_iterate.multipliers.constraints);
      initial_iterate.feasibility_multipliers.reset();

      // reformulate (scale, add slacks, relax the bounds, ...) if necessary
      std::unique_ptr<Model> model = ModelFactory::reformulate(std::move(ampl_model), initial_iterate, options);

      // create the constraint relaxation strategy, the globalization mechanism and the Uno solver
      auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(*model, options);
      auto globalization_mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
      Uno uno = Uno(*globalization_mechanism, options);

      // solve the instance
      Result result = uno.solve(*model, initial_iterate, options);
      Uno::print_optimization_summary(options, result);
      // std::cout << "memory_allocation_amount = " << memory_allocation_amount << '\n';
   }

   Level Logger::level = INFO;
} // namespace

int main(int argc, char* argv[]) {
   using namespace uno;
   
   if (1 < argc) {
      // get the default options
      Options options = Options::get_default_options("uno.options");
      // override them with the command line arguments
      options.get_command_line_arguments(argc, argv);
      Logger::set_logger(options.get_string("logger"));

      if (std::string(argv[1]) == "-v") {
         Uno::print_uno_version();
      }
      else if (std::string(argv[1]) == "--strategies") {
         Uno::print_available_strategies();
      }
      else {
         // run Uno on the .nl file (last command line argument)
         std::string model_name = std::string(argv[argc - 1]);
         run_uno_ampl(model_name, options);
      }
   }
   else {
      Uno::print_uno_version();
   }
   return EXIT_SUCCESS;
}
