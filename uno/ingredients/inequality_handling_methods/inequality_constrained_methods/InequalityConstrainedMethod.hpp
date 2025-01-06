// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include "../InequalityHandlingMethod.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class InequalityConstrainedMethod : public InequalityHandlingMethod {
   public:
      InequalityConstrainedMethod(const std::string& hessian_model, size_t number_variables, size_t number_hessian_nonzeros, bool convexify,
            const Options& options);
      ~InequalityConstrainedMethod() override = default;
      
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void set_initial_point(const Vector<double>& point) override;
      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      [[nodiscard]] double proximal_coefficient(const Iterate& current_iterate) const override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;

      void set_auxiliary_measure(const Model& model, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate&, const Vector<double>&, double) const override;

      void postprocess_iterate(const OptimizationProblem& model, Iterate& iterate) override;

   protected:
      Vector<double> initial_point{};

      static void compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H
