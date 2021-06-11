#ifndef ITERATE_H
#define ITERATE_H

#include <ostream>
#include <vector>
#include "Problem.hpp"
#include "Constraint.hpp"
#include "Matrix.hpp"

enum TerminationStatus {
   NOT_OPTIMAL = 0, KKT_POINT, /* feasible stationary point */
   FJ_POINT, /* infeasible stationary point */
   FEASIBLE_SMALL_STEP, INFEASIBLE_SMALL_STEP
};

struct Residuals {
   double constraints;
   double KKT;
   double FJ;
   double complementarity;
};

struct ProgressMeasures {
   double feasibility;
   double objective;
};

/*! \class Iterate
 * \brief Optimization iterate
 *
 *  Point and its evaluations during an optimization process
 */
class Iterate {
public:
   /*!
    *  Constructor
    */
   Iterate() = default;
   Iterate(size_t number_variables, size_t number_constraints);
   Iterate(const std::vector<double>& x, const Multipliers& multipliers);

   std::vector<double> x; /*!< \f$\mathbb{R}^n\f$ primal variables */
   Multipliers multipliers; /*!< \f$\mathbb{R}^n\f$ Lagrange multipliers/dual variables */
   static int number_eval_objective;
   static int number_eval_constraints;
   static int number_eval_jacobian;

   // functions
   double objective; /*!< Objective value */
   bool is_objective_computed;

   std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
   bool are_constraints_computed;

   SparseVector objective_gradient; /*!< Sparse Jacobian of the objective */
   bool is_objective_gradient_computed; /*!< Flag that indicates if the objective gradient has already been computed */

   std::vector<SparseVector> constraints_jacobian; /*!< Sparse Jacobian of the constraints */
   bool is_constraints_jacobian_computed; /*!< Flag that indicates if the constraint Jacobian has already been computed */

   // residuals
   Residuals residuals;
   //double feasibility_measure;
   //double optimality_measure;
   ProgressMeasures progress;

   void compute_objective(const Problem& problem);
   void compute_constraints(const Problem& problem);
   void compute_objective_gradient(const Problem& problem);
   void set_objective_gradient(const SparseVector& objective_gradient);
   void compute_constraints_jacobian(const Problem& problem);
   std::vector<double> lagrangian_gradient(const Problem& problem, double objective_multiplier, const Multipliers& multipliers);

   void clear();

   friend std::ostream& operator<<(std::ostream& stream, const Iterate& iterate);
};

#endif // ITERATE_H
