#ifndef UNO_QPSUBPROBLEM_H
#define UNO_QPSUBPROBLEM_H

#include "ActiveSetSubproblem.hpp"
#include "HessianModel.hpp"
#include "solvers/QP/QPSolver.hpp"
#include "tools/Options.hpp"

class QPSubproblem : public ActiveSetSubproblem {
public:
   QPSubproblem(const NonlinearReformulation& problem, const Options& options);

   [[nodiscard]] Direction solve(Statistics& statistics, const NonlinearReformulation& problem, Iterate& current_iterate) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const NonlinearReformulation& problem, const Direction& direction) const override;
   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] double get_proximal_coefficient() const override;

protected:
   // use pointers to allow polymorphism
   const std::unique_ptr<QPSolver> solver; /*!< Solver that solves the subproblem */
   const double proximal_coefficient;

   // evaluations
   const std::unique_ptr<HessianModel> hessian_model; /*!< Strategy to evaluate or approximate the Hessian */

   void evaluate_problem(const NonlinearReformulation& problem, Iterate& current_iterate);
};

#endif // UNO_QPSUBPROBLEM_H
