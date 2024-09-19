#include <pybind11/pybind11.h>
#include "linear_algebra/RectangularMatrix.hpp"

namespace py = pybind11;

namespace uno {
   void define_RectangularMatrix(py::module& module) {
      // RectangularMatrix class (instantiated with double elements)
      py::class_<RectangularMatrix<double>>(module, "RectangularMatrix")
         // constructor
         .def(py::init<size_t, size_t>(), py::arg("number_rows"), py::arg("number_columns"), "Constructor")
         // methods
         .def("get_number_rows", &RectangularMatrix<double>::get_number_rows, "Number of rows")
         .def("get_number_columns", &RectangularMatrix<double>::get_number_columns, "Number of columns")
         .def("clear", &RectangularMatrix<double>::clear, "Empty the matrix")
         .def("__getitem__", [](RectangularMatrix<double>& matrix, size_t row_index) -> SparseVector<double>& {
            return matrix[row_index];
         }, py::arg("row_index"), "Returns a row", py::return_value_policy::reference_internal);
   }
} // namespace
