#include <pybind11/pybind11.h>

#include <PhilipsPixelEngine/softwarerendercontext.hpp>

namespace py = pybind11;

PYBIND11_MODULE(softwarerendercontext, m)
{
    py::module_::import("pixelengine");

    py::class_<SoftwareRenderContext, RenderContext>(m, "SoftwareRenderContext")
        .def(py::init<>())
        .def(py::init<size_t, size_t>(), py::arg("width"), py::arg("height"))
        .def_property("number_of_worker_threads",
            py::overload_cast<>(&SoftwareRenderContext::numberOfWorkerThreads, py::const_),
            py::overload_cast<size_t>(&SoftwareRenderContext::numberOfWorkerThreads));
}