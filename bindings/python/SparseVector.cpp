#include <pybind11/pybind11.h>
#include "linear_algebra/SparseVector.hpp"

namespace py = pybind11;

namespace uno {
   void define_SparseVector(py::module& module) {
      // SparseVector class (instantiated with double elements)
      py::class_<SparseVector<double>>(module, "SparseVector")
         // constructor
         .def(py::init<size_t>(), py::arg("capacity"), "Constructor")
         // methods
         .def("size", &SparseVector<double>::size, "Number of elements of the sparse vector")
         .def("insert", &SparseVector<double>::insert, py::arg("index"), py::arg("value"), "Insert an element into the sparse vector")
         .def("clear", &SparseVector<double>::clear, "Delete all elements in the sparse vector")
         .def("is_empty", &SparseVector<double>::is_empty, "True if the sparse vector is empty, false otherwise")
         // string representation
         .def("__repr__", [](const SparseVector<double>& vector) {
            std::string representation = "Sparse vector:";
            for (const auto [index, value]: vector) {
               representation += "\n(" + std::to_string(index) + ", " + std::to_string(value) + ")";
            }
            return representation;
         });
   }
} // namespace
