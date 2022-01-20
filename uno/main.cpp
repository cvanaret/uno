#include "interfaces/AMPL/AMPLModel.hpp"
#include "ingredients/strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategyFactory.hpp"
#include "Uno.hpp"
#include "optimization/ScaledReformulation.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"

// new() overload to track heap allocations
size_t total_allocations = 0;

/*
void* operator new(size_t size) {
   std::cout << "Allocating " << size << " bytes\n";
   total_allocations += size;
   return malloc(size);
}
*/

void run_uno_ampl(const std::string& problem_name, const Options& options) {
   // TODO: use a factory
   // AMPL model
   auto original_problem = std::make_unique<AMPLModel>(problem_name);
   INFO << "Heap allocations after AMPL: " << total_allocations << "\n";

   // initial primal and dual points
   Iterate first_iterate(original_problem->number_variables, original_problem->number_constraints);
   original_problem->get_initial_primal_point(first_iterate.x);
   original_problem->get_initial_dual_point(first_iterate.multipliers.constraints);
   // project x into the bounds
   original_problem->project_point_in_bounds(first_iterate.x);

   // initialize the function scaling TODO put this constant in option file
   Scaling scaling(original_problem->number_constraints, 100.);
   // function scaling
   const bool scale_functions = (options.at("scale_functions") == "yes");
   if (scale_functions) {
      // evaluate the gradients at the current point
      first_iterate.evaluate_objective_gradient(*original_problem);
      first_iterate.evaluate_constraint_jacobian(*original_problem);
      scaling.compute(first_iterate.objective_gradient, first_iterate.constraint_jacobian);
      // forget about these evaluations
      first_iterate.reset_evaluations();
   }
   const Problem& problem_to_solve = ScaledReformulation(*original_problem, scaling);

   // create the constraint relaxation strategy
   auto constraint_relaxation_strategy = ConstraintRelaxationStrategyFactory::create(problem_to_solve, options);
   INFO << "Heap allocations after ConstraintRelax, Subproblem and Solver: " << total_allocations << "\n";

   // create the globalization mechanism
   auto mechanism = GlobalizationMechanismFactory::create(*constraint_relaxation_strategy, options);
   INFO << "Heap allocations after Mechanism: " << total_allocations << "\n";

   Uno uno = Uno(*mechanism, options);

   INFO << "Heap allocations before solving: " << total_allocations << "\n";
   const bool enforce_linear_constraints = (options.at("enforce_linear_constraints") == "yes");
   Result result = uno.solve(problem_to_solve, first_iterate, enforce_linear_constraints);
   Uno::postsolve_solution(*original_problem, scaling, result.solution, result.status);

   const bool print_solution = (options.at("print_solution") == "yes");
   result.print(print_solution);
   INFO << "Heap allocations: " << total_allocations << "\n";
}

Level Logger::logger_level = INFO;

int main(int argc, char* argv[]) {
   if (1 < argc) {
      // get the default options
      Options options = get_default_options("uno.options");
      // get the command line options
      get_command_line_options(argc, argv, options);
      set_logger(options.at("logger"));

      options.print();

      if (std::string(argv[1]) == "-v") {
         std::cout << "Welcome in Uno\n";
         std::cout << "To solve an AMPL problem, type ./uno_ampl path_to_file/file.nl\n";
         std::cout << "To choose a globalization mechanism, use the argument -mechanism [LS|TR]\n";
         std::cout << "To choose a globalization strategy, use the argument -strategy [penalty|filter|nonmonotone-filter]\n";
         std::cout << "To choose a constraint relaxation strategy, use the argument -constraint-relaxation [feasibility-restoration|l1-relaxation]\n";
         std::cout << "To choose a subproblem, use the argument -subproblem [QP|LP|barrier]\n";
         std::cout << "To choose a preset, use the argument -preset [byrd|filtersqp|ipopt]\n";
         std::cout << "The options can be combined in the same command line. Autocompletion is active.\n";
      }
      else {
         // run Uno on the .nl file (last command line argument)
         std::string problem_name = std::string(argv[argc - 1]);
         run_uno_ampl(problem_name, options);
      }
   }
   return EXIT_SUCCESS;
}
