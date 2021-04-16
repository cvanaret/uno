#include "Problem.hpp"
#include <cmath>
#include <iostream>
#include "Utils.hpp"

/* Abstract Problem class */

Problem::Problem(std::string name, int number_variables, int number_constraints):
name(name), number_variables(number_variables), number_constraints(number_constraints),
// allocate all vectors
variable_name(number_variables),
//variable_discrete(number_variables),
variables_bounds(number_variables), variable_status(number_variables),
constraint_name(number_constraints),
//constraint_variables(number_constraints),
constraint_bounds(number_constraints), constraint_type(number_constraints), constraint_status(number_constraints),
hessian_maximum_number_nonzeros(0) {
}

Problem::~Problem() {
}

/* compute ||c|| */
double Problem::compute_constraint_residual(std::vector<double>& constraints, std::string norm_value) {
    std::vector<double> residuals(constraints.size());
    for (int j = 0; j < this->number_constraints; j++) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, norm_value);
}

/* compute ||c_S|| for a given set S */
double Problem::compute_constraint_residual(std::vector<double>& constraints, std::set<int> constraint_set, std::string norm_value) {
    SparseGradient residuals;
    for (int j: constraint_set) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, norm_value);
}

void Problem::determine_bounds_types(std::vector<Range>& bounds, std::vector<ConstraintType>& status) {
    for (unsigned int i = 0; i < bounds.size(); i++) {
        if (bounds[i].lb == bounds[i].ub) {
            status[i] = EQUAL_BOUNDS;
        }
        else if (-INFINITY < bounds[i].lb && bounds[i].ub < INFINITY) {
            status[i] = BOUNDED_BOTH_SIDES;
        }
        else if (-INFINITY < bounds[i].lb) {
            status[i] = BOUNDED_LOWER;
        }
        else if (bounds[i].ub < INFINITY) {
            status[i] = BOUNDED_UPPER;
        }
        else {
            status[i] = UNBOUNDED;
        }
    }
    return;
}

void Problem::determine_constraints_() {
    int current_equality_constraint = 0;
    int current_inequality_constraint = 0;
    for (int j = 0; j < this->number_constraints; j++) {
        if (this->constraint_status[j] == EQUAL_BOUNDS) {
            this->equality_constraints[j] = current_equality_constraint;
            current_equality_constraint++;
        }
        else {
            this->inequality_constraints[j] = current_inequality_constraint;
            current_inequality_constraint++;
        }
    }
    return;
}

/* native C++ problem */

//CppProblem::CppProblem(std::string name, int number_variables, int number_constraints, double (*objective)(std::vector<double> x), std::vector<double> (*objective_gradient)(std::vector<double> x)):
//Problem(name, number_variables, number_constraints),
//objective_(objective),
//objective_gradient_(objective_gradient) {
//}
//
//double CppProblem::objective(std::vector<double>& x) {
//    return this->objective_(x);
//}
//
//std::vector<double> CppProblem::objective_dense_gradient(std::vector<double>& x) {
//    return this->objective_gradient_(x);
//}
//
//SparseGradient CppProblem::objective_sparse_gradient(std::vector<double>& x) {
//    std::vector<double> dense_gradient = this->objective_gradient_(x);
//    SparseGradient sparse_gradient;
//    for (unsigned int i = 0; i < dense_gradient.size(); i++) {
//        if (dense_gradient[i] != 0.) {
//            sparse_gradient[i] = dense_gradient[i];
//        }
//    }
//    return sparse_gradient;
//}
//
//double CppProblem::evaluate_constraint(int j, std::vector<double>& x) {
//    return this->constraints_[j](x);
//}
//
//std::vector<double> CppProblem::evaluate_constraints(std::vector<double>& x) {
//    std::vector<double> constraints(this->number_constraints);
//    for (int j = 0; j < this->number_constraints; j++) {
//        constraints[j] = this->evaluate_constraint(j, x);
//    }
//    return constraints;
//}
