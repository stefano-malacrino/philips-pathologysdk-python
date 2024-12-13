#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PhilipsPixelEngine/pixelengine.hpp>
#include <PhilipsPixelEngine/renderbackend.hpp>
#include <PhilipsPixelEngine/rendercontext.hpp>

#include "pystreambuf.h"

namespace py = pybind11;

std::ios_base::openmode parse_isyntax_open_options(std::string const &open_mode, std::string const &container_name) {
    if (!(container_name == "" || container_name == "ficom" || container_name == "caching-ficom"))
            throw py::value_error("invalid container: " + container_name);
    if (open_mode == "r")
        return std::ios::in | std::ios::binary;
    else if (open_mode == "w")
        return std::ios::out | std::ios::binary;
    else
        throw py::value_error("mode needs to be either 'r' or 'w'");
}

PYBIND11_MODULE(pixelengine, m)
{
    py::class_<RenderBackend> pyRenderBackend(m, "RenderBackend");

    py::enum_<RenderBackend::ImageFormatType>(pyRenderBackend, "ImageFormatType")
        .value("RGB", RenderBackend::ImageFormatType::RGB)
        .value("RGBA", RenderBackend::ImageFormatType::RGBA)
        .value("LUMINANCE", RenderBackend::ImageFormatType::LUMINANCE)
        .value("UNDEFINED_FORMAT", RenderBackend::ImageFormatType::UNDEFINED_FORMAT);

    py::class_<RenderContext>(m, "RenderContext")
        .def(py::init<>())
        .def_property_readonly("width", &RenderContext::width)
        .def_property_readonly("height", &RenderContext::height);

    py::class_<PixelEngine::FilterHandle>(m, "FilterHandle");
        //.def("supported_parameters", &PixelEngine::FilterHandle::supportedParameters);

    py::class_<PixelEngine> pyPixelEngine(m, "PixelEngine");

    py::enum_<PixelEngine::BufferType>(pyPixelEngine, "BufferType")
        .value("RGB", PixelEngine::BufferType::RGB)
        .value("RGBA", PixelEngine::BufferType::RGBA)
        .value("LUMINANCE", PixelEngine::BufferType::LUMINANCE);

    py::class_<PixelEngine::CompressionParameters>(pyPixelEngine, "CompressionParameters");

    py::class_<PixelEngine::WSICompressionParametersBuilder>(pyPixelEngine, "WSICompressionParametersBuilder")
        .def("with_num_derived_levels", &PixelEngine::WSICompressionParametersBuilder::withNumDerivedLevels)
        .def("with_block_size", &PixelEngine::WSICompressionParametersBuilder::withBlockSize)
        .def("with_pixel_transform", &PixelEngine::WSICompressionParametersBuilder::withPixelTransform)
        .def("with_compressor", &PixelEngine::WSICompressionParametersBuilder::withCompressor)
        .def("with_colorspace_transform", &PixelEngine::WSICompressionParametersBuilder::withColorspaceTransform)
        .def("with_bit_depth", &PixelEngine::WSICompressionParametersBuilder::withBitDepth)
        .def("with_num_threads", &PixelEngine::WSICompressionParametersBuilder::withNumThreads)
        .def("with_default_color_intensity", &PixelEngine::WSICompressionParametersBuilder::withDefaultColorIntensity)
        //.def("with_quality_preset", &PixelEngine::WSICompressionParametersBuilder::withQualityPreset)
        .def("with_quality", &PixelEngine::WSICompressionParametersBuilder::withQuality)
        .def("with_origin", &PixelEngine::WSICompressionParametersBuilder::withOrigin)
        .def("with_scale", &PixelEngine::WSICompressionParametersBuilder::withScale)
        .def("build", &PixelEngine::WSICompressionParametersBuilder::build);

    py::class_<PixelEngine::SecondaryCaptureCompressionParametersBuilder>(pyPixelEngine, "SecondaryCaptureCompressionParametersBuilder")
        .def("with_origin", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::withOrigin)
        .def("with_scale", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::withScale)
        .def("build", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::build);

    py::class_<PixelEngine::DataEnvelopes>(pyPixelEngine, "DataEnvelopes")
        .def("__eq__", [](PixelEngine::DataEnvelopes &self, PixelEngine::DataEnvelopes &other)
             { return self == other; })
        .def("as_extreme_vertices_model", &PixelEngine::DataEnvelopes::asEVM)
        .def("as_rectangles", &PixelEngine::DataEnvelopes::asRectangles);

    py::class_<PixelEngine::Region, std::shared_ptr<PixelEngine::Region>>(pyPixelEngine, "Region")
        .def("__hash__", &PixelEngine::Region::id)
        .def("__eq__", [](PixelEngine::Region &self, PixelEngine::Region &other)
             { return self.id() == other.id(); })
        .def("__ge__", [](PixelEngine::Region &self, PixelEngine::Region &other)
             { return self.id() >= other.id(); })
        .def("__gt__", [](PixelEngine::Region &self, PixelEngine::Region &other)
             { return self.id() > other.id(); })
        .def("__le__", [](PixelEngine::Region &self, PixelEngine::Region &other)
             { return self.id() <= other.id(); })
        .def("__lt__", [](PixelEngine::Region &self, PixelEngine::Region &other)
             { return self.id() < other.id(); })
        .def_property_readonly("ready", &PixelEngine::Region::ready)
        .def_property_readonly("range", &PixelEngine::Region::range)
        .def("get", [](PixelEngine::Region &self, py::buffer buffer)
             {
            py::buffer_info info = buffer.request();
            return self.get(info.ptr, info.itemsize * info.size); })
        .def("draw", &PixelEngine::Region::draw);

    py::class_<PixelEngine::Level>(pyPixelEngine, "Level")
        .def("chain_source_view", &PixelEngine::Level::chainSourceView)
        .def("filter_width", &PixelEngine::Level::filterWidth);

    py::class_<PixelEngine::View>(pyPixelEngine, "View")
        .def(
            "__getitem__", [](PixelEngine::View &self, size_t level)
            { return &self[level]; },
            py::return_value_policy::reference)
        .def("chain_source_view", [](PixelEngine::View &self, const PixelEngine::View &view, int x_shift, int y_shift, int level_shift)
             { return self.chainSourceView(view, x_shift, y_shift, level_shift); })
        .def(
            "request_regions", [](PixelEngine::View &self, std::vector<std::vector<size_t>> const &regions, bool enable_async_rendering, std::vector<size_t> const &background_color, PixelEngine::BufferType buffer_type)
            { return self.requestRegions(regions, enable_async_rendering, background_color, buffer_type); },
            py::arg("regions"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def(
            "request_regions", [](PixelEngine::View &self, std::vector<std::vector<size_t>> const &regions, PixelEngine::DataEnvelopes const &data_envelopes, bool enable_async_rendering, std::vector<size_t> const &background_color, PixelEngine::BufferType buffer_type)
            { return self.requestRegions(regions, data_envelopes, enable_async_rendering, background_color, buffer_type); },
            py::arg("regions"),
            py::arg("data_envelopes"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def("dimension_ranges", [](PixelEngine::View &self, size_t level)
             { return self.dimensionRanges(level); })
        .def_property_readonly("dimension_names", [](PixelEngine::View &self)
                               { return self.dimensionNames(); })
        .def_property_readonly("dimension_units", [](PixelEngine::View &self)
                               { return self.dimensionUnits(); })
        .def_property_readonly("dimension_types", [](PixelEngine::View &self)
                               { return self.dimensionTypes(); })
        .def_property_readonly("dimension_discrete_values", [](PixelEngine::View &self)
                               { return self.dimensionDiscreteValues(); })
        .def_property_readonly("scale", [](PixelEngine::View &self)
                               { return self.scale(); })
        .def_property_readonly("origin", [](PixelEngine::View &self)
                               { return self.origin(); })
        .def(
            "data_envelopes", [](PixelEngine::View &self, size_t level)
            { return &self.dataEnvelopes(level); },
            py::return_value_policy::reference)
        .def_property_readonly("bits_allocated", [](PixelEngine::View &self)
                               { return self.bitsAllocated(); })
        .def_property_readonly("bits_stored", [](PixelEngine::View &self)
                               { return self.bitsStored(); })
        .def_property_readonly("high_bit", [](PixelEngine::View &self)
                               { return self.highBit(); })
        .def_property_readonly("pixel_representation", [](PixelEngine::View &self)
                               { return self.pixelRepresentation(); })
        .def_property_readonly("planar_configuration", [](PixelEngine::View &self)
                               { return self.planarConfiguration(); })
        .def_property_readonly("samples_per_pixel", [](PixelEngine::View &self)
                               { return self.samplesPerPixel(); })
        .def_property_readonly("id", [](PixelEngine::View &self)
                               { return self.id(); })
        .def_property_readonly("num_derived_levels", [](PixelEngine::View &self)
                               { return self.numDerivedLevels(); })
        .def_property_readonly("pixel_size", [](PixelEngine::View &self)
                               { return self.pixelSize(); });

    py::class_<PixelEngine::SourceView>(pyPixelEngine, "SourceView")
        .def(
            "__getitem__", [](PixelEngine::SourceView &self, size_t level)
            { return &self[level]; },
            py::return_value_policy::reference)
        .def("chain_source_view", [](PixelEngine::SourceView &self, const PixelEngine::View &view, int x_shift, int y_shift, int level_shift)
             { return self.chainSourceView(view, x_shift, y_shift, level_shift); })
        .def(
            "request_regions", [](PixelEngine::SourceView &self, std::vector<std::vector<size_t>> const &regions, bool enable_async_rendering, std::vector<size_t> const &background_color, PixelEngine::BufferType buffer_type)
            { return self.requestRegions(regions, enable_async_rendering, background_color, buffer_type); },
            py::arg("regions"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def(
            "request_regions", [](PixelEngine::SourceView &self, std::vector<std::vector<size_t>> const &regions, PixelEngine::DataEnvelopes const &data_envelopes, bool enable_async_rendering, std::vector<size_t> const &background_color, PixelEngine::BufferType buffer_type)
            { return self.requestRegions(regions, data_envelopes, enable_async_rendering, background_color, buffer_type); },
            py::arg("regions"),
            py::arg("data_envelopes"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def("dimension_ranges", [](PixelEngine::SourceView &self, size_t level)
             { return self.dimensionRanges(level); })
        .def_property_readonly("dimension_names", [](PixelEngine::SourceView &self)
                               { return self.dimensionNames(); })
        .def_property_readonly("dimension_units", [](PixelEngine::SourceView &self)
                               { return self.dimensionUnits(); })
        .def_property_readonly("dimension_types", [](PixelEngine::SourceView &self)
                               { return self.dimensionTypes(); })
        .def_property_readonly("dimension_discrete_values", [](PixelEngine::SourceView &self)
                               { return self.dimensionDiscreteValues(); })
        .def_property_readonly("scale", [](PixelEngine::SourceView &self)
                               { return self.scale(); })
        .def_property_readonly("origin", [](PixelEngine::SourceView &self)
                               { return self.origin(); })
        .def(
            "data_envelopes", [](PixelEngine::SourceView &self, size_t level)
            { return &self.dataEnvelopes(level); },
            py::return_value_policy::reference)
        .def_property_readonly("bits_allocated", [](PixelEngine::SourceView &self)
                               { return self.bitsAllocated(); })
        .def_property_readonly("bits_stored", [](PixelEngine::SourceView &self)
                               { return self.bitsStored(); })
        .def_property_readonly("high_bit", [](PixelEngine::SourceView &self)
                               { return self.highBit(); })
        .def_property_readonly("pixel_representation", [](PixelEngine::SourceView &self)
                               { return self.pixelRepresentation(); })
        .def_property_readonly("planar_configuration", [](PixelEngine::SourceView &self)
                               { return self.planarConfiguration(); })
        .def_property_readonly("samples_per_pixel", [](PixelEngine::SourceView &self)
                               { return self.samplesPerPixel(); })
        .def_property_readonly("id", [](PixelEngine::SourceView &self)
                               { return self.id(); })
        .def_property_readonly("num_derived_levels", [](PixelEngine::SourceView &self)
                               { return self.numDerivedLevels(); })
        .def_property_readonly("pixel_size", [](PixelEngine::SourceView &self)
                               { return self.pixelSize(); })
        .def("load_default_parameters", &PixelEngine::SourceView::loadDefaultParameters)
        .def("truncation", &PixelEngine::SourceView::truncation);

    py::class_<PixelEngine::DisplayView, PixelEngine::SourceView>(pyPixelEngine, "DisplayView")
        .def_property(
            "sharpness", [](PixelEngine::DisplayView &self)
            { return self.sharpness(); },
            [](PixelEngine::DisplayView &self, double gain)
            { return self.sharpness(gain); })
        .def_property(
            "contrast_clip_limit", [](PixelEngine::DisplayView &self)
            { return self.contrastClipLimit(); },
            [](PixelEngine::DisplayView &self, double clip_limit)
            { return self.contrastClipLimit(clip_limit); })
        .def_property(
            "color_correction_gamma", [](PixelEngine::DisplayView &self)
            { return self.colorCorrectionGamma(); },
            [](PixelEngine::DisplayView &self, double gamma)
            { return self.colorCorrectionGamma(gamma); })
        .def_property(
            "color_correction_black_point", [](PixelEngine::DisplayView &self)
            { return self.colorCorrectionBlackPoint(); },
            [](PixelEngine::DisplayView &self, double blackpoint)
            { return self.colorCorrectionBlackPoint(blackpoint); })
        .def_property(
            "color_correction_white_point", [](PixelEngine::DisplayView &self)
            { return self.colorCorrectionWhitePoint(); },
            [](PixelEngine::DisplayView &self, double whitepoint)
            { return self.colorCorrectionWhitePoint(whitepoint); })
        .def_property(
            "color_gain", [](PixelEngine::DisplayView &self)
            { return self.colorGain(); },
            [](PixelEngine::DisplayView &self, double gain)
            { return self.colorGain(gain); });

    py::class_<PixelEngine::SubImage>(pyPixelEngine, "SubImage")
        .def("include_input_region", &PixelEngine::SubImage::includeInputRegion)
        .def("preallocate_pixels", &PixelEngine::SubImage::preallocatePixels)
        /*.def("put_block", [](PixelEngine::SubImage &self, py::buffer buffer)
             { 
                py::buffer_info info = buffer.request();
                return self.putBlock(info.ptr, info.itemsize * info.size); })*/
        .def_property_readonly("has_display_view", &PixelEngine::SubImage::hasDisplayView)
        .def_property_readonly(
            "source_view", [](PixelEngine::SubImage &self)
            { return &self.sourceView(); },
            py::return_value_policy::reference)
        .def_property_readonly(
            "display_view", [](PixelEngine::SubImage &self)
            { return &self.displayView(); },
            py::return_value_policy::reference)
        .def("add_view", &PixelEngine::SubImage::addView)
        .def(
            "block_size", [](PixelEngine::SubImage &self, size_t template_id)
            { return self.blockSize(template_id); },
            py::arg("template_id") = 0UL)
        .def("block_pos", &PixelEngine::SubImage::blockPos)
        .def("read_block", [](PixelEngine::SubImage &self, py::buffer buffer)
             { 
                py::buffer_info info = buffer.request();
                return self.readBlock(info.ptr, info.itemsize * info.size); })
        //.def("ordered_block_coordinates", &PixelEngine::SubImage::orderedBlockCoordinates)
        .def_property_readonly("image_type", &PixelEngine::SubImage::imageType)
        .def_property_readonly("pixel_transform", &PixelEngine::SubImage::pixelTransform)
        //.def_property_readonly("quality_preset", &PixelEngine::SubImage::qualityPreset)
        .def_property_readonly("quality", &PixelEngine::SubImage::quality)
        .def_property_readonly("compressor", &PixelEngine::SubImage::compressor)
        .def_property_readonly("colorspace_transform", &PixelEngine::SubImage::colorspaceTransform)
        .def_property_readonly("num_tiles", &PixelEngine::SubImage::numTiles)
        .def_property(
            "icc_profile", [](PixelEngine::SubImage &self)
            { return self.iccProfile(); },
            [](PixelEngine::SubImage &self, std::string const &profile)
            { return self.iccProfile(profile); })
        //.def_property_readonly("icc_matrix", &PixelEngine::SubImage::iccMatrix)
        .def_property(
            "image_data", [](PixelEngine::SubImage &self)
            {
        std::vector<uint8_t> const &data = self.imageData();
        return py::memoryview::from_memory((void*) data.data(), data.size(), true);},
            [](PixelEngine::SubImage &self, py::buffer const &data)
            {
        py::buffer_info info = data.request();
        std::vector<uint8_t> values((uint8_t *)info.ptr, ((uint8_t *)info.ptr) + info.itemsize * info.size);
        return self.imageData(values); })
        .def_property(
            "lossy_image_compression", [](PixelEngine::SubImage &self)
            { return self.lossyImageCompression(); },
            [](PixelEngine::SubImage &self, std::string const &val)
            { return self.lossyImageCompression(val); })
        .def_property(
            "lossy_image_compression_ratio", [](PixelEngine::SubImage &self)
            { return self.lossyImageCompressionRatio(); },
            [](PixelEngine::SubImage &self, double val)
            { return self.lossyImageCompressionRatio(val); })
        .def_property(
            "lossy_image_compression_method", [](PixelEngine::SubImage &self)
            { return self.lossyImageCompressionMethod(); },
            [](PixelEngine::SubImage &self, std::string const &method)
            { return self.lossyImageCompressionMethod(method); })
        .def_property(
            "color_linearity", [](PixelEngine::SubImage &self)
            { return self.colorLinearity(); },
            [](PixelEngine::SubImage &self, std::string const &color_linearity)
            { return self.colorLinearity(color_linearity); })
        //.def("header", &PixelEngine::SubImage::header)
        //.def("block_coordinate", &PixelEngine::SubImage::blockCoordinate)
        //.def("block_index", &PixelEngine::SubImage::blockIndex)
        ;

    py::class_<PixelEngine::ISyntaxFacade>(pyPixelEngine, "ISyntaxFacade")
        .def(
            "open", [](PixelEngine::ISyntaxFacade &self, std::string const &url, std::string const &container_name, std::string const &mode, std::string const &cache_name)
            {
        std::ios_base::openmode open_mode = parse_isyntax_open_options(mode, container_name);
        return self.open(url, container_name, open_mode, cache_name); },
            py::arg("url"),
            py::arg("container_name") = "",
            py::arg("mode") = "r",
            py::arg("cache_name") = "")
        .def(
            "open", [](PixelEngine::ISyntaxFacade &self, std::istream& stream, std::string const &container_name, std::string mode, std::string const &cache_name)
            {
        std::ios_base::openmode open_mode = parse_isyntax_open_options(mode, container_name);
        return self.open(static_cast<std::iostream*>(&stream), container_name, open_mode, cache_name); },
            py::arg("stream"),
            py::arg("container_name") = "",
            py::arg("mode") = "r",
            py::arg("cache_name") = "")
        .def("add_sub_image", &PixelEngine::ISyntaxFacade::addSubImage)
        .def_property_readonly("num_images", &PixelEngine::ISyntaxFacade::numImages)
        .def(
            "__getitem__", [](PixelEngine::ISyntaxFacade &self, size_t index)
            { return &self[index]; },
            py::return_value_policy::reference)
        .def(
            "__getitem__", [](PixelEngine::ISyntaxFacade &self, std::string const &type)
            { return &self[type]; },
            py::return_value_policy::reference)
        .def("finalize_geometry_and_properties", &PixelEngine::ISyntaxFacade::finalizeGeometryAndProperties)
        .def_property_readonly("isyntax_file_version", &PixelEngine::ISyntaxFacade::iSyntaxFileVersion)
        .def_property_readonly("id", &PixelEngine::ISyntaxFacade::id)
        .def_property(
            "barcode", [](PixelEngine::ISyntaxFacade &self)
            { return self.barcode(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &code)
            { return self.barcode(code); })
        .def_property(
            "scanner_calibration_status", [](PixelEngine::ISyntaxFacade &self)
            { return self.scannerCalibrationStatus(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &status)
            { return self.scannerCalibrationStatus(status); })
        .def_property(
            "software_versions", [](PixelEngine::ISyntaxFacade &self)
            { return self.softwareVersions(); },
            [](PixelEngine::ISyntaxFacade &self, std::vector<std::string> const &versions)
            { return self.softwareVersions(versions); })
        .def_property(
            "derivation_description", [](PixelEngine::ISyntaxFacade &self)
            { return self.derivationDescription(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &description)
            { return self.derivationDescription(description); })
        .def_property(
            "acquisition_datetime", [](PixelEngine::ISyntaxFacade &self)
            { return self.acquisitionDateTime(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &datetime)
            { return self.acquisitionDateTime(datetime); })
        .def_property(
            "manufacturer", [](PixelEngine::ISyntaxFacade &self)
            { return self.manufacturer(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &make)
            { return self.manufacturer(make); })
        .def_property(
            "model_name", [](PixelEngine::ISyntaxFacade &self)
            { return self.modelName(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &model)
            { return self.modelName(model); })
        .def_property(
            "device_serial_number", [](PixelEngine::ISyntaxFacade &self)
            { return self.deviceSerialNumber(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &serial)
            { return self.deviceSerialNumber(serial); })
        .def_property(
            "scanner_rack_number", [](PixelEngine::ISyntaxFacade &self)
            { return self.scannerRackNumber(); },
            [](PixelEngine::ISyntaxFacade &self, uint16_t rack_number)
            { return self.scannerRackNumber(rack_number); })
        .def_property(
            "scanner_slot_number", [](PixelEngine::ISyntaxFacade &self)
            { return self.scannerSlotNumber(); },
            [](PixelEngine::ISyntaxFacade &self, uint16_t slot_number)
            { return self.scannerSlotNumber(slot_number); })
        .def_property(
            "scanner_operator_id", [](PixelEngine::ISyntaxFacade &self)
            { return self.scannerOperatorId(); },
            [](PixelEngine::ISyntaxFacade &self, std::string const &operator_id)
            { return self.scannerOperatorId(operator_id); })
        .def_property(
            "scanner_rack_priority", [](PixelEngine::ISyntaxFacade &self)
            { return self.scannerRackPriority(); },
            [](PixelEngine::ISyntaxFacade &self, uint16_t priority)
            { return self.scannerRackPriority(priority); })
        .def_property(
            "date_of_last_calibration", [](PixelEngine::ISyntaxFacade &self)
            { return self.dateOfLastCalibration(); },
            [](PixelEngine::ISyntaxFacade &self, std::vector<std::string> const &date)
            { return self.dateOfLastCalibration(date); })
        .def_property(
            "time_of_last_calibration", [](PixelEngine::ISyntaxFacade &self)
            { return self.timeOfLastCalibration(); },
            [](PixelEngine::ISyntaxFacade &self, std::vector<std::string> const &time)
            { return self.timeOfLastCalibration(time); })
        .def_property_readonly("is_philips", &PixelEngine::ISyntaxFacade::isPhilips)
        //.def_property_readonly("is_hamamatsu", &PixelEngine::ISyntaxFacade::isHamamatsu)
        .def_property_readonly("is_UFS", &PixelEngine::ISyntaxFacade::isUFS)
        .def_property_readonly("is_UFSb", &PixelEngine::ISyntaxFacade::isUFSb)
        .def_property_readonly("is_UVS", &PixelEngine::ISyntaxFacade::isUVS)
        .def("close", &PixelEngine::ISyntaxFacade::close)
        .def("abort", &PixelEngine::ISyntaxFacade::abort)
        //.def(".remaining_pixels_to_encode", &PixelEngine::ISyntaxFacade::remainingPixelsToEncode)
        ;

    pyPixelEngine.def(py::init<>())
        .def(py::init<RenderBackend &, RenderContext &>())
        .def_property_readonly_static("version", [](py::object)
                                      { return PixelEngine::version(); })
        .def(
            "__getitem__", [](PixelEngine &self, std::string const &name)
            { return &self[name]; },
            py::return_value_policy::reference)
        .def("containers", &PixelEngine::containers)
        .def("container_version", &PixelEngine::containerVersion)
        .def("compressors", &PixelEngine::compressors)
        .def("pixel_transforms", &PixelEngine::pixelTransforms)
        .def("block_sizes", &PixelEngine::blockSizes)
        .def("colorspace_transforms", &PixelEngine::colorspaceTransforms)
        .def("quality_presets", &PixelEngine::qualityPresets)
        //.def("supported_filters", &PixelEngine::supportedFilters)
        .def("wait_all", &PixelEngine::waitAll)
        .def(
            "wait_any", [](PixelEngine &self)
            { return self.waitAny(); },
            py::return_value_policy::reference)
        .def(
            "wait_any", [](PixelEngine &self, std::list<std::shared_ptr<PixelEngine::Region>> const &regions)
            { return self.waitAny(regions); },
            py::return_value_policy::reference)
        .def("clear_render_target", &PixelEngine::clearRenderTarget)
        .def("clear_render_cache", &PixelEngine::clearRenderCache)
        .def("clear_render_buffers", &PixelEngine::clearRenderBuffers)
        .def_property(
            "network_timeout", [](PixelEngine &self)
            { return self.networkTimeout(); },
            [](PixelEngine &self, size_t timeout)
            { return self.networkTimeout(timeout); })
        .def(
            "client_certificates", [](PixelEngine &self, std::string const &cert, std::string const &key, std::string const &password)
            { return self.clientCertificates(cert, key, password); },
            py::arg("cert"),
            py::arg("key") = "",
            py::arg("password") = "")
        .def_property(
            "certificates", [](PixelEngine &self)
            { return self.certificates(); },
            [](PixelEngine &self, std::string const &path)
            { return self.certificates(path); })
        .def_property(
            "rate_limiting", [](PixelEngine &self)
            { return self.rateLimiting(); },
            [](PixelEngine &self, bool enabled)
            { return self.rateLimiting(enabled); });
}