// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSEVECTOR_H
#define UNO_SPARSEVECTOR_H

#include <cassert>
#include <functional>
#include "tools/Logger.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Collection.hpp"

// SparseVector is a sparse vector that uses contiguous memory. It contains:
// - a vector of indices of type size_t
// - a vector of values of type ElementType
// the indices are neither unique nor sorted
template <typename ElementType>
class SparseVector: public Collection<ElementType> {
public:
   class iterator {
   public:
      iterator(const SparseVector<ElementType>& vector, size_t index) : vector(vector), index(index) {
      }
      std::pair<size_t, ElementType> operator*() {
         // copy the element in the pair. Cheap only for trivial types
         return {this->vector.indices[this->index], this->vector.values[this->index]};
      }
      // prefix increment
      iterator& operator++() { this->index++; return *this; }

      friend bool operator==(const iterator& a, const iterator& b) { return &(a.vector) == &(b.vector) && a.index == b.index; };
      friend bool operator!=(const iterator& a, const iterator& b) { return &(a.vector) != &(b.vector) || a.index != b.index; };

   private:
      const SparseVector<ElementType>& vector;
      size_t index;
   };

   explicit SparseVector(size_t capacity = 0);
   void for_each(const std::function<void(size_t, ElementType)>& f) const override;
   // void for_each_index(const std::function<void (size_t)>& f) const;
   void for_each_value(const std::function<void(ElementType)>& f) const;

   [[nodiscard]] iterator begin() const { return iterator(*this, 0); }
   [[nodiscard]] iterator end() const { return iterator(*this, this->number_nonzeros); }

   [[nodiscard]] size_t size() const;
   void reserve(size_t capacity);

   void insert(size_t index, ElementType value);
   void transform(const std::function<ElementType(ElementType)>& f);
   void clear();
   [[nodiscard]] bool empty() const;

   template <typename U>
   friend std::ostream& operator<<(std::ostream& stream, const SparseVector<U>& x);

protected:
   std::vector<size_t> indices{};
   std::vector<ElementType> values{};
   size_t number_nonzeros{0};
};

// SparseVector methods
template <typename ElementType>
SparseVector<ElementType>::SparseVector(size_t capacity): Collection<ElementType>() {
   this->reserve(capacity);
}

template <typename ElementType>
void SparseVector<ElementType>::for_each(const std::function<void(size_t, ElementType)>& f) const {
   for (size_t index: Range(this->number_nonzeros)) {
      f(this->indices[index], this->values[index]);
   }
}

template <typename ElementType>
void SparseVector<ElementType>::for_each_value(const std::function<void(ElementType)>& f) const {
   for (size_t index: Range(this->number_nonzeros)) {
      f(this->values[index]);
   }
}

template <typename ElementType>
size_t SparseVector<ElementType>::size() const {
   return this->number_nonzeros;
}

template <typename ElementType>
void SparseVector<ElementType>::reserve(size_t capacity) {
   this->indices.reserve(capacity);
   this->values.reserve(capacity);
}

template <typename ElementType>
void SparseVector<ElementType>::insert(size_t index, ElementType value) {
   this->indices.push_back(index);
   this->values.push_back(value);
   this->number_nonzeros++;
}

template <typename ElementType>
void SparseVector<ElementType>::clear() {
   this->indices.clear();
   this->values.clear();
   this->number_nonzeros = 0;
}

template <typename ElementType>
bool SparseVector<ElementType>::empty() const {
   return (this->number_nonzeros == 0);
}

template <typename ElementType>
void SparseVector<ElementType>::transform(const std::function<ElementType (ElementType)>& f) {
   for (size_t index: Range(this->number_nonzeros)) {
      this->values[index] = f(this->values[index]);
   }
}

template <typename ElementType>
std::ostream& operator<<(std::ostream& stream, const SparseVector<ElementType>& x) {
   stream << "sparse vector with " << x.size() << " nonzeros\n";
   for (const auto [index, element]: x) {
      stream << "index " << index << ", value " << element << '\n';
   }
   return stream;
}

// free functions

template <typename ElementType>
ElementType norm_inf(const SparseVector<ElementType>& x) {
   ElementType norm = ElementType(0);
   for (const auto [index, element]: x) {
      norm = std::max(norm, std::abs(element));
   }
   return norm;
}

template <typename ElementType>
ElementType dot(const std::vector<ElementType>& x, const SparseVector<ElementType>& y) {
   ElementType dot_product = ElementType(0);
   for (const auto [index, y_element]: y) {
      assert(index < x.size() && "Vector.dot: the sparse vector y is larger than the dense vector x");
      dot_product += x[index] * y_element;
   }
   return dot_product;
}

// precondition: factor != 0
template <typename ElementType>
void scale(SparseVector<ElementType>& x, ElementType factor) {
   x.transform([=](ElementType element) {
      return factor * element;
   });
}

#endif // UNO_SPARSEVECTOR_H
