#ifndef LPSOLVERFACTORY_H
#define LPSOLVERFACTORY_H

#include <memory>
#include "LPSolver.hpp"

#ifdef HAS_BQPD
#include "BQPDSolver.hpp"
#endif

class LPSolverFactory {
public:
   static std::unique_ptr<BQPDSolver> create(size_t number_variables, size_t number_constraints, const std::string& LP_solver_name) {
#ifdef HAS_BQPD
      if (LP_solver_name == "BQPD") {
         return std::make_unique<BQPDSolver>(number_variables, number_constraints, 0, false);
      }
      throw std::invalid_argument("LP solver not found");
#endif
   }
};

#endif // LPSOLVERFACTORY_H
