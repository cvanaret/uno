#include <pybind11/pybind11.h>
#include "linear_algebra/SparseVector.hpp"

namespace py = pybind11;

namespace uno {
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python binding to the solver Uno for nonconvex optimization";
      
      // SparseVector class (instantiated with double elements)
      py::class_<SparseVector<double>>(module, "SparseVector")
         // constructor
         .def(py::init<size_t>(), py::arg("capacity"))
         // methods
         .def("size", &SparseVector<double>::size, "Number of elements of the sparse vector")
         .def("insert", &SparseVector<double>::insert, "Insert an element into the sparse vector")
         .def("clear", &SparseVector<double>::clear, "Delete all elements in the sparse vector")
         .def("is_empty", &SparseVector<double>::is_empty, "True if the sparse vector is empty, false otherwise")
         // string representation
         .def("__repr__", [](const SparseVector<double>& vector) {
            std::string representation = "Sparse vector:\n";
            for (const auto [index, value]: vector) {
               representation += "(" + std::to_string(index) + ", " + std::to_string(value) + ")\n";
            }
            return representation;
         })
      ;
   }
} // namespace
