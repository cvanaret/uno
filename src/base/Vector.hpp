#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
#include "Logger.hpp"
#include "SparseGradient.hpp"

enum Norm {L1_NORM = 1, L2_NORM = 2, L2_SQUARED_NORM, INF_NORM};

std::vector<double> add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor = 1.);

double norm_1(const std::vector<double>& x);
double norm_1(const SparseGradient& x);
double norm_1(const std::vector<SparseGradient>& m);

double norm_2_squared(const std::vector<double>& x);
double norm_2_squared(const SparseGradient& x);
double norm_2(const std::vector<double>& x);
double norm_2(const SparseGradient& x);

double norm_inf(const std::vector<double>& x, size_t length = std::numeric_limits<unsigned int>::max());
double norm_inf(const SparseGradient& x);
double norm_inf(const std::vector<SparseGradient>& m);

double dot(const std::vector<double>& x, const std::vector<double>& y);
double dot(const std::vector<double>& x, const SparseGradient& y);
double dot(const SparseGradient& x, const SparseGradient& y);

template <typename T>
double norm(const T& x, Norm norm) {
    /* choose the right norm */
    if (norm == INF_NORM) {
        return norm_inf(x);
    }
    else if (norm == L2_NORM) {
        return norm_2(x);
    }
    else if (norm == L2_SQUARED_NORM) {
        return norm_2_squared(x);
    }
    else if (norm == L1_NORM) {
        return norm_1(x);
    }
    else {
        throw std::out_of_range("The norm is not known");
    }
}

template <typename T>
void print_vector(std::ostream &stream, const std::vector<T>& x, const char end='\n', unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (size_t i = start; i < std::min<unsigned int>(start + length, x.size()); i++) {
        stream << x[i] << " ";
    }
    stream << end;
    return;
}

template <typename T>
void print_vector(const Level& level, const std::vector<T>& x, const char end='\n', unsigned int start = 0, unsigned int length = std::numeric_limits<unsigned int>::max()) {
    for (size_t i = start; i < std::min<unsigned int>(start + length, x.size()); i++) {
        level << x[i] << " ";
    }
    level << end;
    return;
}

template <typename T>
void print_vector(const Level& level, const std::set<T>& x, const char end='\n') {
    for (T xi: x) {
        level << xi << " ";
    }
    level << end;
    return;
}

template <typename T, typename U>
void print_vector(std::ostream &stream, const std::map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        stream << "x[" << i << "] = " << xi << ", ";
    }
    stream << end;
    return;
}

template <typename T, typename U>
void print_vector(const Level& level, const std::map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        level << "x[" << i << "] = " << xi << ", ";
    }
    level << end;
    return;
}

template <typename T, typename U>
void print_vector(std::ostream &stream, const std::unordered_map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        stream << "x[" << i << "] = " << xi << ", ";
    }
    stream << end;
    return;
}

template <typename T, typename U>
void print_vector(const Level& level, const std::unordered_map<T, U>& x, const char end='\n') {
    for (std::pair<T, U> element : x) {
        T i = element.first;
        U xi = element.second;
        level << "x[" << i << "] = " << xi << ", ";
    }
    level << end;
    return;
}

std::string join(std::vector<std::string>& vector, const std::string& separator);

#endif // UTILS_H
