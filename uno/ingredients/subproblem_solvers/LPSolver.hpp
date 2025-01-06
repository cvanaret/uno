// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LPSOLVER_H
#define UNO_LPSOLVER_H

namespace uno {
   // forward declarations
   class Direction;
   class LagrangeNewtonSubproblem;
   class Statistics;
   template <typename ElementType>
   class Vector;
   struct WarmstartInformation;

   /*! \class LPSolver
    * \brief LP solver
    *
    */
   class LPSolver {
   public:
      LPSolver() = default;
      virtual ~LPSolver() = default;

      virtual void solve_LP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) = 0;
   };
} // namespace

#endif // UNO_LPSOLVER_H
