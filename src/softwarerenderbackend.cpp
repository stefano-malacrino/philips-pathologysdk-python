#include <pybind11/pybind11.h>

#include <PhilipsPixelEngine/softwarerenderbackend.hpp>

namespace py = pybind11;

PYBIND11_MODULE(softwarerenderbackend, m)
{
    py::module_::import("pixelengine");
    
    py::class_<SoftwareRenderBackend, RenderBackend>(m, "SoftwareRenderBackend")
        .def(py::init<RenderBackend::ImageFormatType, size_t>(),
            py::arg("image_format_type") = RenderBackend::ImageFormatType::RGBA,
            py::arg("tile_cache_size_bytes") = 256000000);
}