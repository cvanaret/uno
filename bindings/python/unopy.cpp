#include <pybind11/pybind11.h>
#include "linear_algebra/SparseVector.hpp"

namespace py = pybind11;

namespace uno {
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python binding to the solver Uno for nonconvex optimization";
      
      // SparseVector class (instantiated with double elements)
      py::class_<SparseVector<double>>(module, "SparseVector")
         .def(py::init<size_t>())
         .def("size", &SparseVector<double>::size, "")
         .def("insert", &SparseVector<double>::insert, "")
         .def("clear", &SparseVector<double>::clear, "")
         .def("is_empty", &SparseVector<double>::is_empty, "")
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
