#include "QPSolverFactory.hpp"
#ifdef BQPD
#include "BQPDSolver.hpp"
#endif

std::shared_ptr<QPSolver> QPSolverFactory::create(const std::string& QP_solver_name, int number_variables, int number_constraints, int maximum_number_nonzeros, bool quadratic_programming) {
#ifdef BQPD
    if (QP_solver_name == "BQPD") {
        return std::make_shared<BQPDSolver>(number_variables, number_constraints, maximum_number_nonzeros, quadratic_programming);
    }
#endif
    throw std::invalid_argument("QPSolver name " + QP_solver_name + " does not exist");
}
