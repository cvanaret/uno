// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOUNDRELAXEDMODEL_H
#define UNO_BOUNDRELAXEDMODEL_H

#include <memory>
#include "Model.hpp"

namespace uno {
   // forward declaration
   class Options;

   class BoundRelaxedModel: public Model {
   public:
      BoundRelaxedModel(std::unique_ptr<Model> original_model, const Options& options);

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override { return this->model->evaluate_objective(x); }
      void evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const override {
         this->model->evaluate_objective_gradient(x, gradient);
      }
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override {
         this->model->evaluate_constraints(x, constraints);
      }
      void evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override {
         this->model->evaluate_constraint_gradient(x, constraint_index, gradient);
      }
      void evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override {
         this->model->evaluate_constraint_jacobian(x, constraint_jacobian);
      }
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            SymmetricMatrix<size_t, double>& hessian) const override {
         this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
      }
      void compute_hessian_vector_product(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            Vector<double>& result) const override {
         this->model->compute_hessian_vector_product(x, objective_multiplier, multipliers, result);
      }

      // only these two functions are redefined
      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override { return this->model->get_variable_bound_type(variable_index); }
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->model->get_lower_bounded_variables(); }
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->model->get_upper_bounded_variables(); }
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->model->get_slacks(); }
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->model->get_single_lower_bounded_variables(); }
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->model->get_single_upper_bounded_variables(); }
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override { return this->model->get_fixed_variables(); }

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override { return this->model->constraint_lower_bound(constraint_index); }
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override { return this->model->constraint_upper_bound(constraint_index); }
      [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override { return this->model->get_constraint_type(constraint_index); }
      [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override { return this->model->get_constraint_bound_type(constraint_index); }
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->model->get_equality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->model->get_inequality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->model->get_linear_constraints(); }

      void initial_primal_point(Vector<double>& x) const override { this->model->initial_primal_point(x); }
      void initial_dual_point(Vector<double>& multipliers) const override { this->model->initial_dual_point(multipliers); }
      void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const override {
         this->model->postprocess_solution(iterate, termination_status);
      }

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override { return this->model->number_objective_gradient_nonzeros(); }
      [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model->number_jacobian_nonzeros(); }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model->number_hessian_nonzeros(); }

   private:
      const std::unique_ptr<Model> model{};
      const double relaxation_factor;
   };
} // namespace

#endif // UNO_BOUNDRELAXEDMODEL_H