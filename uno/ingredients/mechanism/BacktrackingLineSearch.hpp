#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "GlobalizationMechanism.hpp"

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class BacktrackingLineSearch : public GlobalizationMechanism {
public:
   explicit BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations = 7, double backtracking_ratio
   = 0.5);

   void initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) override;
   std::tuple<Iterate, double> compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate& current_iterate) override;

private:
   double step_length{1.};
   /* ratio of step length update in ]0, 1[ */
   const double backtracking_ratio;
   const double min_step_length{1e-6};

   bool termination_();
   void print_iteration_();
   void add_statistics(Statistics& statistics, const Direction& direction);
   void decrease_step_length();
};

#endif // LINESEARCH_H
