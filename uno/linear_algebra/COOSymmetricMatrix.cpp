#include <iostream>
#include "COOSymmetricMatrix.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOSymmetricMatrix::COOSymmetricMatrix(size_t dimension, size_t capacity) : SymmetricMatrix(dimension, capacity) {
   row_indices.reserve(capacity);
   column_indices.reserve(capacity);
}

void COOSymmetricMatrix::reset() {
   SymmetricMatrix::reset();
   this->matrix.clear();
   this->row_indices.clear();
   this->column_indices.clear();
}

// generic iterator
void COOSymmetricMatrix::for_each(const std::function<void(size_t, size_t, double)>& f) const {
   for (size_t k = 0; k < this->number_nonzeros; k++) {
      const size_t i = this->row_indices[k];
      const size_t j = this->column_indices[k];
      const double entry = this->matrix[k];
      f(i, j, entry);
   }
}

size_t COOSymmetricMatrix::find(size_t row_index, size_t column_index) const {
   size_t position = 0;
   while (position < this->number_nonzeros) {
      if (this->row_indices[position] == row_index && this->column_indices[position] == column_index) {
         break;
      }
      position++;
   }
   return position;
}

void COOSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   // look for an existing entry
   size_t position = this->find(row_index, column_index);
   // if the entry already exists, modify it
   if (false && position < this->number_nonzeros) {
      this->matrix[position] += term;
   }
   // otherwise, add an entry
   else {
      this->matrix.push_back(term);
      this->row_indices.push_back(row_index);
      this->column_indices.push_back(column_index);
      this->number_nonzeros++;
   }
}

void COOSymmetricMatrix::pop() {
   this->matrix.pop_back();
   this->row_indices.pop_back();
   this->column_indices.pop_back();
   this->number_nonzeros--;
}

void COOSymmetricMatrix::add_identity_multiple(double multiple) {
   for (size_t i = 0; i < this->dimension; i++) {
      this->insert(multiple, i, i);
   }
}

double COOSymmetricMatrix::smallest_diagonal_entry() const {
   double smallest_entry = std::numeric_limits<double>::infinity();
   this->for_each([&](size_t i, size_t j, double entry) {
      if (i == j) {
         smallest_entry = std::min(smallest_entry, entry);
      }
   });
   if (smallest_entry == std::numeric_limits<double>::infinity()) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

void COOSymmetricMatrix::print(std::ostream& stream) const {
   this->for_each([&](size_t i, size_t j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << "\n";
   });
}