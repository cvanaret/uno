// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_QPSOLVER_H
#define UNO_QPSOLVER_H

#include <vector>
#include "LPSolver.hpp"

namespace uno {
   // forward declarations
   class Direction;
   //class HessianModel;
   class LagrangeNewtonSubproblem;
   class Statistics;
   struct WarmstartInformation;

   class QPSolver : public LPSolver {
   public:
      QPSolver(): LPSolver() { }
      ~QPSolver() override = default;

      void solve_LP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override = 0;

      virtual void solve_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) = 0;

      [[nodiscard]] virtual double hessian_quadratic_product(const Vector<double>& primal_direction) const = 0;
   };
} // namespace

#endif // UNO_QPSOLVER_H