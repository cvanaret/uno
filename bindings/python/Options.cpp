#include <pybind11/pybind11.h>
#include "tools/Options.hpp"

namespace py = pybind11;

namespace uno {
   void define_Options(py::module& module) {
      // Options class
      py::class_<Options>(module, "Options")
         // constructor
         .def(py::init<>(), "Constructor")
         // methods
         .def("__setitem__", [](Options& options, const std::string& key, const std::string& value) {
            options[key] = value;
         }, py::arg("key"), py::arg("value"), "Set an option")
         .def("__getitem__", [](const Options& options, const std::string& key) {
            return options.at(key);
         }, py::arg("key"), "Read an option")
         .def_static("get_default_options", &Options::get_default_options, py::arg("file_path"), "Read default options from an option file")
         // string representation
         .def("__repr__", [](const Options& options) {
            return options.to_string(false);
         })
      ;
   }
} // namespace