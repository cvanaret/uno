#ifndef HESSIANEVALUATION_H
#define HESSIANEVALUATION_H

#include <memory>
#include <vector>
#include "optimization/Problem.hpp"
#include "solvers/linear/LinearSolver.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"

// virtual (abstract) class
template<class SYMMETRIC_MATRIX>
class HessianModel {
public:
   explicit HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros);
   virtual ~HessianModel() = default;

   SYMMETRIC_MATRIX hessian;
   int evaluation_count{0};

   virtual void evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) = 0;

   SYMMETRIC_MATRIX modify_inertia(SYMMETRIC_MATRIX& matrix, LinearSolver<SYMMETRIC_MATRIX>& linear_solver);
};

template<class SYMMETRIC_MATRIX>
inline HessianModel<SYMMETRIC_MATRIX>::HessianModel(size_t dimension, size_t hessian_maximum_number_nonzeros):
      hessian(dimension, hessian_maximum_number_nonzeros) {
}

template<class SYMMETRIC_MATRIX>
inline SYMMETRIC_MATRIX HessianModel<SYMMETRIC_MATRIX>::modify_inertia(SYMMETRIC_MATRIX& matrix, LinearSolver<SYMMETRIC_MATRIX>& linear_solver) {
   const double beta = 1e-4;

   // Nocedal and Wright, p51
   double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the Hessian is " << matrix.smallest_diagonal_entry() << "\n";

   double inertia = 0.;
   double previous_inertia = 0.;
   if (smallest_diagonal_entry <= 0.) {
      inertia = beta - smallest_diagonal_entry;
   }

   if (0. < inertia) {
      matrix.add_identity_multiple(inertia - previous_inertia);
   }

   DEBUG << "Testing factorization with inertia term " << inertia << "\n";
   linear_solver.do_symbolic_factorization(matrix);
   linear_solver.do_numerical_factorization(matrix);

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << linear_solver.number_negative_eigenvalues() << " negative eigenvalues\n";
      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == 0) {
         good_inertia = true;
         DEBUG << "Factorization was a success with inertia " << inertia << "\n";
      }
      else {
         previous_inertia = inertia;
         inertia = (inertia == 0.) ? beta : 2 * inertia;
         matrix.add_identity_multiple(inertia - previous_inertia);
         DEBUG << "Testing factorization with inertia term " << inertia << "\n";
         linear_solver.do_numerical_factorization(matrix);
      }
   }
   return matrix;
}

// Exact Hessian

template<class SYMMETRIC_MATRIX>
class ExactHessian : public HessianModel<SYMMETRIC_MATRIX> {
public:
   explicit ExactHessian(size_t dimension, size_t hessian_maximum_number_nonzeros) : HessianModel<SYMMETRIC_MATRIX>(dimension,
         hessian_maximum_number_nonzeros) {
   }

   ~ExactHessian() override = default;

   void evaluate(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
         const std::vector<double>& constraint_multipliers) override {
      /* compute Hessian */
      problem.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
      this->evaluation_count++;
   }
};

//template <class SYMMETRIC_MATRIX>
//class ConvexifiedHessian : public HessianEvaluation<SYMMETRIC_MATRIX> {
//public:
//   ConvexifiedHessian(size_t dimension, size_t hessian_maximum_number_nonzeros, const std::string& linear_solver_name):
//      HessianEvaluation<SYMMETRIC_MATRIX>(dimension, hessian_maximum_number_nonzeros), linear_solver(LinearSolverFactory<MA57Solver>::create
//      (linear_solver_name)) {
//   }
//   ~ConvexifiedHessian() override = default;
//
//   void compute(const Problem& problem, const std::vector<double>& primal_variables, double objective_multiplier,
//         const std::vector<double>& constraint_multipliers) override {
//      /* compute Hessian */
//      problem.lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, this->hessian);
//      this->evaluation_count++;
//      /* modify the inertia to make the problem strictly convex */
//      DEBUG << "hessian before convexification: " << this->hessian;
//      this->hessian = this->modify_inertia(this->hessian, *this->linear_solver);
//   }
//
//protected:
//   std::unique_ptr<LinearSolver<COOSymmetricMatrix>> linear_solver; /*!< Solver that computes the inertia */
//};

template<class SYMMETRIC_MATRIX>
class HessianModelFactory {
public:
   static std::unique_ptr <HessianModel<SYMMETRIC_MATRIX>> create(const std::string& hessian_model_strategy, size_t dimension,
         size_t hessian_maximum_number_nonzeros, bool convexify) {
      if (hessian_model_strategy == "exact") {
         if (convexify) {
            //   return std::make_unique<ConvexifiedHessian<SYMMETRIC_MATRIX>>(dimension, hessian_maximum_number_nonzeros, "MA57");
         }
         //else {
         return std::make_unique<ExactHessian<SYMMETRIC_MATRIX>>(dimension, hessian_maximum_number_nonzeros);
         //}
      }
      else {
         throw std::invalid_argument("Hessian evaluation method " + hessian_model_strategy + " does not exist");
      }
   }
};

#endif // HESSIANEVALUATION_H