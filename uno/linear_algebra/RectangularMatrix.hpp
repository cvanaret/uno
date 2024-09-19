// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RECTANGULARMATRIX_H
#define UNO_RECTANGULARMATRIX_H

#include <vector>
#include "SparseVector.hpp"

namespace uno {
   // TODO use more appropriate sparse representation
   //template <typename ElementType>
   //using RectangularMatrix = std::vector<SparseVector<ElementType>>;

   template <typename ElementType>
   class RectangularMatrix {
   public:
      using value_type = ElementType;

      RectangularMatrix(size_t number_rows, size_t number_columns): matrix(number_rows), number_rows(number_rows), number_columns(number_columns) {
         for (auto& row: this->matrix) {
            row.reserve(number_columns);
         }
      }

      [[nodiscard]] size_t get_number_rows() const {
         return this->number_rows;
      }

      [[nodiscard]] size_t get_number_columns() const {
         return this->number_columns;
      }

      SparseVector<ElementType>& operator[](size_t row_index) {
         return this->matrix[row_index];
      }

      const SparseVector<ElementType>& operator[](size_t row_index) const {
         return this->matrix[row_index];
      }

      void clear() {
         for (auto& row: this->matrix) {
            row.clear();
         }
      }

   protected:
      std::vector<SparseVector<ElementType>> matrix;
      size_t number_rows;
      size_t number_columns;
   };
} // namespace

#endif // UNO_RECTANGULARMATRIX_H
