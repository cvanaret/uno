// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDSOLVER_H
#define UNO_BQPDSOLVER_H

#include <array>
#include <vector>
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "ingredients/subproblem_solvers/QPSolver.hpp"

namespace uno {
   // forward declarations
   class Multipliers;
   class Options;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

   // see bqpd.f
   enum class BQPDStatus {
      OPTIMAL = 0,
      UNBOUNDED_PROBLEM = 1,
      BOUND_INCONSISTENCY = 2,
      INFEASIBLE = 3,
      INCORRECT_PARAMETER = 4,
      LP_INSUFFICIENT_SPACE = 5,
      HESSIAN_INSUFFICIENT_SPACE = 6,
      SPARSE_INSUFFICIENT_SPACE = 7,
      MAX_RESTARTS_REACHED = 8,
      UNDEFINED = 9
   };

   enum BQPDMode {
      COLD_START = 0,
      ACTIVE_SET_EQUALITIES = 1, // cold start
      USER_DEFINED = 2, // hot start
      UNCHANGED_ACTIVE_SET = 3,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN = 4,
      UNCHANGED_ACTIVE_SET_AND_REDUCED_HESSIAN = 5,
      UNCHANGED_ACTIVE_SET_AND_JACOBIAN_AND_REDUCED_HESSIAN = 6, // warm start
   };

   enum class BQPDProblemType {LP, QP};

   class BQPDSolver : public QPSolver {
   public:
      BQPDSolver(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, BQPDProblemType problem_type, const Options& options);

      void solve_LP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

      void solve_QP(Statistics& statistics, LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& primal_direction) const override;

   private:
      std::vector<double> lower_bounds{}, upper_bounds{}; // lower and upper bounds of variables and constraints
      std::vector<double> constraints;
      SparseVector<double> linear_objective;
      RectangularMatrix<double> constraint_jacobian;
      std::vector<double> bqpd_jacobian{};
      std::vector<int> bqpd_jacobian_sparsity{};
      SymmetricMatrix<size_t, double> hessian;

      int kmax{0}, mlp{1000};
      size_t mxwk0{2000000}, mxiwk0{500000};
      std::array<int, 100> info{};
      std::vector<double> alp{};
      std::vector<int> lp{}, active_set{};
      std::vector<double> w{}, gradient_solution{}, residuals{}, e{};
      size_t size_hessian_sparsity{};
      size_t size_hessian_workspace{};
      size_t size_hessian_sparsity_workspace{};
      std::vector<double> workspace{};
      std::vector<int> workspace_sparsity{};
      int k{0};
      int iprint{0}, nout{6};
      double fmin{-1e20};
      int peq_solution{0}, ifail{0};
      const int fortran_shift{1};
      Vector<int> current_hessian_indices{};

      const bool print_subproblem;

      void set_up_subproblem(LagrangeNewtonSubproblem& subproblem, const WarmstartInformation& warmstart_information);
      void solve_subproblem(LagrangeNewtonSubproblem& subproblem, const Vector<double>& initial_point, Direction& direction,
            const WarmstartInformation& warmstart_information);
      [[nodiscard]] static BQPDMode determine_mode(const WarmstartInformation& warmstart_information);
      void save_hessian_to_local_format();
      void save_gradients_to_local_format(size_t number_constraints);
      void set_multipliers(size_t number_variables, Multipliers& direction_multipliers);
      static BQPDStatus bqpd_status_from_int(int ifail);
      static SubproblemStatus status_from_bqpd_status(BQPDStatus bqpd_status);
   };
} // namespace

#endif // UNO_BQPDSOLVER_H