// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFEASIBLEINTERIORPOINTMETHOD_H
#define UNO_INFEASIBLEINTERIORPOINTMETHOD_H

#include "../InequalityHandlingMethod.hpp"
#include "PrimalDualInteriorPointProblem.hpp"
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "BarrierParameterUpdateStrategy.hpp"

namespace uno {
   // forward references
   template <typename IndexType, typename NumericalType>
   class DirectSymmetricIndefiniteLinearSolver;
   class DualResiduals;

   struct InteriorPointParameters {
      double tau_min;
      double k_sigma;
      double regularization_exponent;
      double small_direction_factor;
      double push_variable_to_interior_k1;
      double push_variable_to_interior_k2;
   };

   class PrimalDualInteriorPointMethod : public InequalityHandlingMethod {
   public:
      PrimalDualInteriorPointMethod(size_t number_variables, size_t number_constraints, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, const Options& options);

      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void set_initial_point(const Vector<double>& point) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& constraint_index) override;
      [[nodiscard]] double proximal_coefficient(const Iterate& current_iterate) const override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;

      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
            Direction& direction, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& primal_direction) const override;

      void set_auxiliary_measure(const Model& model, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate& current_iterate,
            const Vector<double>& primal_direction, double step_length) const override;

      void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) override;

   protected:
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      Vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      RectangularMatrix<double> constraint_jacobian; /*!< Sparse Jacobian of the constraints */
      SymmetricMatrix<size_t, double> hessian;

      SymmetricIndefiniteLinearSystem<double> augmented_system;
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver;

      BarrierParameterUpdateStrategy barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const double default_multiplier;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double damping_factor; // (Section 3.7 in IPOPT paper)
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool solving_feasibility_problem{false};
      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      [[nodiscard]] double push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const;
      void evaluate_functions(Statistics& statistics, const PrimalDualInteriorPointProblem& barrier_problem, Iterate& current_iterate,
            const Multipliers& current_multipliers, const WarmstartInformation& warmstart_information);
      void update_barrier_parameter(const OptimizationProblem& problem, const Iterate& current_iterate, const Multipliers& current_multipliers,
            const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
      [[nodiscard]] double compute_barrier_term_directional_derivative(const Model& model, const Iterate& current_iterate,
            const Vector<double>& primal_direction) const;
      [[nodiscard]] static double primal_fraction_to_boundary(const OptimizationProblem& problem, const Vector<double>& current_primals,
            const Vector<double>& primal_direction, double tau);
      [[nodiscard]] static double dual_fraction_to_boundary(const OptimizationProblem& problem, const Multipliers& current_multipliers,
            Multipliers& direction_multipliers, double tau);
      void assemble_augmented_system(Statistics& statistics, const OptimizationProblem& problem, const Multipliers& current_multipliers,
            WarmstartInformation& warmstart_information);
      void assemble_augmented_rhs(const Multipliers& current_multipliers, size_t number_variables, size_t number_constraints);
      void assemble_primal_dual_direction(const OptimizationProblem& problem, const Vector<double>& current_primals, const Multipliers& current_multipliers,
            Vector<double>& direction_primals, Multipliers& direction_multipliers);
      void compute_bound_dual_direction(const OptimizationProblem& problem, const Vector<double>& current_primals, const Multipliers& current_multipliers,
            const Vector<double>& primal_direction, Multipliers& direction_multipliers);
      void compute_least_square_multipliers(const OptimizationProblem& problem, Iterate& iterate, Vector<double>& constraint_multipliers);
   };
} // namespace

#endif // UNO_INFEASIBLEINTERIORPOINTMETHOD_H