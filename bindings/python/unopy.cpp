#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace uno {
   // individual classes
   void define_SparseVector(py::module& module);
   void define_RectangularMatrix(py::module& module);
   void define_Options(py::module& module);

   // module definition
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python binding to the solver Uno for nonconvex optimization";
      
      define_SparseVector(module);
      define_RectangularMatrix(module);
      define_Options(module);
   }
} // namespace
