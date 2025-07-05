#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <PhilipsPixelEngine/pixelengine.hpp>
#include <PhilipsPixelEngine/renderbackend.hpp>
#include <PhilipsPixelEngine/rendercontext.hpp>

#include "pystreambuf.h"

namespace py = pybind11;

std::ios_base::openmode parse_isyntax_open_mode(std::string const &open_mode) {
    if (open_mode == "r")
        return std::ios::in | std::ios::binary;
    else if (open_mode == "w")
        return std::ios::out | std::ios::binary;
    else
        throw py::value_error("mode needs to be either 'r' or 'w'");
}

class FilterHandle
{
public:
    FilterHandle(std::unique_ptr<PixelEngine::FilterHandle> filter_handle) {
        this->filter_handle = std::move(filter_handle);
    }
    std::unique_ptr<PixelEngine::FilterHandle> &ref() {
        return filter_handle;
    }
    std::string const& name() {
        return filter_handle->name();
    }
    std::vector<std::string> const& supportedParameters() {
        return filter_handle->supportedParameters();
    }

private:
    std::unique_ptr<PixelEngine::FilterHandle> filter_handle;
};

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

    py::class_<PixelEngine> pyPixelEngine(m, "PixelEngine");

    py::enum_<PixelEngine::BufferType>(pyPixelEngine, "BufferType")
        .value("RGB", PixelEngine::BufferType::RGB)
        .value("RGBA", PixelEngine::BufferType::RGBA)
        .value("LUMINANCE", PixelEngine::BufferType::LUMINANCE);

    py::class_<PixelEngine::CompressionParameters, std::shared_ptr<PixelEngine::CompressionParameters>>(pyPixelEngine, "CompressionParameters");

    py::class_<PixelEngine::WSICompressionParametersBuilder>(pyPixelEngine, "WSICompressionParametersBuilder")
        .def(py::init<std::vector<size_t> const&>(), py::arg("size"))
        .def("with_num_derived_levels", &PixelEngine::WSICompressionParametersBuilder::withNumDerivedLevels, py::arg("num_levels"))
        .def("with_block_size", &PixelEngine::WSICompressionParametersBuilder::withBlockSize, py::arg("block_size"))
        .def("with_pixel_transform", &PixelEngine::WSICompressionParametersBuilder::withPixelTransform, py::arg("pixel_transform"))
        .def("with_compressor", &PixelEngine::WSICompressionParametersBuilder::withCompressor, py::arg("compressor"))
        .def("with_colorspace_transform", &PixelEngine::WSICompressionParametersBuilder::withColorspaceTransform, py::arg("colorspace_transform"))
        .def("with_bit_depth", &PixelEngine::WSICompressionParametersBuilder::withBitDepth, py::arg("bit_depth"))
        .def("with_num_threads", &PixelEngine::WSICompressionParametersBuilder::withNumThreads, py::arg("num_threads"))
        .def("with_default_color_intensity", &PixelEngine::WSICompressionParametersBuilder::withDefaultColorIntensity, py::arg("default_color_intensity"))
        .def("with_quality_preset", &PixelEngine::WSICompressionParametersBuilder::withQualityPreset, py::arg("quality"))
        .def("with_quality", &PixelEngine::WSICompressionParametersBuilder::withQuality, py::arg("quality"))
        .def("with_origin", &PixelEngine::WSICompressionParametersBuilder::withOrigin, py::arg("origin"))
        .def("with_scale", &PixelEngine::WSICompressionParametersBuilder::withScale, py::arg("scale"))
        .def("build", &PixelEngine::WSICompressionParametersBuilder::build);
    
    py::class_<PixelEngine::SecondaryCaptureCompressionParametersBuilder>(pyPixelEngine, "SecondaryCaptureCompressionParametersBuilder")
        .def(py::init<std::string const&, uint32_t, uint32_t>(), py::arg("image_type"), py::arg("width"), py::arg("height"))
        .def("with_origin", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::withOrigin, py::arg("origin"))
        .def("with_scale", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::withScale, py::arg("scale"))
        .def("build", &PixelEngine::SecondaryCaptureCompressionParametersBuilder::build);

    py::class_<PixelEngine::DataEnvelopes>(pyPixelEngine, "DataEnvelopes")
        .def("__eq__", &PixelEngine::DataEnvelopes::operator==, py::arg("other"))
        .def("as_extreme_vertices_model", &PixelEngine::DataEnvelopes::asEVM)
        .def("as_rectangles", &PixelEngine::DataEnvelopes::asRectangles);

    py::class_<PixelEngine::Region, std::shared_ptr<PixelEngine::Region>>(pyPixelEngine, "Region")
        .def("__hash__", &PixelEngine::Region::id)
        .def("__eq__", [](PixelEngine::Region &self, PixelEngine::Region const &other)
             { return self.id() == other.id(); }, py::arg("other"))
        .def("__ge__", [](PixelEngine::Region &self, PixelEngine::Region const &other)
             { return self.id() >= other.id(); }, py::arg("other"))
        .def("__gt__", [](PixelEngine::Region &self, PixelEngine::Region const &other)
             { return self.id() > other.id(); }, py::arg("other"))
        .def("__le__", [](PixelEngine::Region &self, PixelEngine::Region const &other)
             { return self.id() <= other.id(); }, py::arg("other"))
        .def("__lt__", [](PixelEngine::Region &self, PixelEngine::Region const &other)
             { return self.id() < other.id(); }, py::arg("other"))
        .def_property_readonly("ready", &PixelEngine::Region::ready)
        .def_property_readonly("range", &PixelEngine::Region::range)
        .def("get", [](PixelEngine::Region &self, py::buffer &buffer)
             {
            py::buffer_info info = buffer.request();
            return self.get(info.ptr, info.itemsize * info.size); },
            py::arg("buffer"))
        .def("draw", &PixelEngine::Region::draw, py::arg("target") = 0);

    py::class_<FilterHandle>(pyPixelEngine, "FilterHandle")
        .def("name", &FilterHandle::name)
        .def("supported_parameters", &FilterHandle::supportedParameters);

    py::class_<PixelEngine::Level>(pyPixelEngine, "Level")
        .def("chain_source_view", &PixelEngine::Level::chainSourceView,
            py::arg("view"), py::arg("x_shift"), py::arg("y_shift"), py::arg("level_shift"))
        .def("filter_width", &PixelEngine::Level::filterWidth,
            py::arg("coord"), py::arg("dimensions"), py::arg("filter_kernel_half_width"));

    py::class_<PixelEngine::View>(pyPixelEngine, "View", py::multiple_inheritance())
        .def("__getitem__",
            &PixelEngine::View::operator[], py::arg("level"),
            py::return_value_policy::reference)
        .def("chain_source_view",
            &PixelEngine::View::chainSourceView,
             py::arg("view"), py::arg("x_shift"), py::arg("y_shift"), py::arg("level_shift"))
        .def("add_user_view", &PixelEngine::View::addChainedView, py::return_value_policy::reference)
        .def("request_regions",
            py::overload_cast<std::vector<std::vector<size_t>> const&, bool, std::vector<size_t> const&, PixelEngine::BufferType>(&PixelEngine::View::requestRegions),
            py::arg("regions"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def("request_regions",
            py::overload_cast<std::vector<std::vector<size_t>> const&, PixelEngine::DataEnvelopes const&, bool, std::vector<size_t> const&, PixelEngine::BufferType>(&PixelEngine::View::requestRegions),
            py::arg("regions"),
            py::arg("data_envelopes"),
            py::arg("enable_async_rendering") = true,
            py::arg("background_color") = std::vector<size_t>({0, 0, 0}),
            py::arg("buffer_type") = PixelEngine::BufferType::RGB,
            py::return_value_policy::reference)
        .def("dimension_ranges", &PixelEngine::View::dimensionRanges, py::arg("level"))
        .def_property_readonly("dimension_names", &PixelEngine::View::dimensionNames)
        .def_property_readonly("dimension_units", &PixelEngine::View::dimensionUnits)
        .def_property_readonly("dimension_types", &PixelEngine::View::dimensionTypes)
        .def_property_readonly("dimension_discrete_values", &PixelEngine::View::dimensionDiscreteValues)
        .def_property_readonly("scale", &PixelEngine::View::scale)
        .def_property_readonly("origin", &PixelEngine::View::origin)
        .def("data_envelopes", &PixelEngine::View::dataEnvelopes, py::arg("level"),
            py::return_value_policy::reference)
        .def_property_readonly("bits_allocated", &PixelEngine::View::bitsAllocated)
        .def_property_readonly("bits_stored", &PixelEngine::View::bitsStored)
        .def_property_readonly("high_bit", &PixelEngine::View::highBit)
        .def_property_readonly("pixel_representation", &PixelEngine::View::pixelRepresentation)
        .def_property_readonly("planar_configuration", &PixelEngine::View::planarConfiguration)
        .def_property_readonly("samples_per_pixel", &PixelEngine::View::samplesPerPixel)
        .def_property_readonly("id", &PixelEngine::View::id)
        .def_property_readonly("num_derived_levels", &PixelEngine::View::numDerivedLevels)
        .def_property_readonly("pixel_size", &PixelEngine::View::pixelSize);

    py::class_<PixelEngine::SourceView>(pyPixelEngine, "SourceView")
        .def(
            "__getitem__", [](PixelEngine::SourceView &self, size_t level)
            { return &self[level]; },
            py::arg("level"),
            py::return_value_policy::reference)
        .def("chain_source_view", [](PixelEngine::SourceView &self, const PixelEngine::View &view, int x_shift, int y_shift, int level_shift)
             { return self.chainSourceView(view, x_shift, y_shift, level_shift); },
            py::arg("view"), py::arg("x_shift"), py::arg("y_shift"), py::arg("level_shift"))
        .def("add_user_view", [](PixelEngine::SourceView &self)
             { return &self.addChainedView(); }, py::return_value_policy::reference)
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
             { return self.dimensionRanges(level); },
             py::arg("level"))
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
            py::arg("level"),
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
        .def("truncation", &PixelEngine::SourceView::truncation,
        py::arg("enabled"), py::arg("rounding"), py::arg("trunc_levels"));

    py::class_<PixelEngine::DisplayView, PixelEngine::SourceView>(pyPixelEngine, "DisplayView")
        .def_property("sharpness",
            py::overload_cast<>(&PixelEngine::DisplayView::sharpness, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::sharpness))
        .def_property("contrast_clip_limit",
            py::overload_cast<>(&PixelEngine::DisplayView::contrastClipLimit, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::contrastClipLimit))
        .def_property("color_correction_gamma",
            py::overload_cast<>(&PixelEngine::DisplayView::colorCorrectionGamma, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::colorCorrectionGamma))
        .def_property("color_correction_black_point",
            py::overload_cast<>(&PixelEngine::DisplayView::colorCorrectionBlackPoint, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::colorCorrectionBlackPoint))
        .def_property("color_correction_white_point",
            py::overload_cast<>(&PixelEngine::DisplayView::colorCorrectionWhitePoint, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::colorCorrectionWhitePoint))
        .def_property("color_gain",
            py::overload_cast<>(&PixelEngine::DisplayView::colorGain, py::const_),
            py::overload_cast<double>(&PixelEngine::DisplayView::colorGain));

    py::class_<PixelEngine::UserView, PixelEngine::View>(pyPixelEngine, "UserView")
        .def("add_filter", [](PixelEngine::UserView &self, std::string const& filter_name){
                return FilterHandle(self.addFilter(filter_name));
            }, py::arg("filter_name"), py::return_value_policy::take_ownership)
        .def("filter_parameter_double", [](PixelEngine::UserView &self, FilterHandle& filter, std::string const& name, double value){
                std::unique_ptr<PixelEngine::FilterHandle> &filter_handle = filter.ref();
                return self.filterParameterDouble(filter_handle, name, value);
            }, py::arg("filter"), py::arg("name"), py::arg("value"))
        .def("filter_parameter_matrix3x3", [](PixelEngine::UserView &self, FilterHandle& filter, std::string const& name, std::array<double, 9> const& value){
                std::unique_ptr<PixelEngine::FilterHandle> &filter_handle = filter.ref();
                return self.filterParameterMatrix3x3(filter_handle, name, value);
            }, py::arg("filter"), py::arg("name"), py::arg("value"));

    py::class_<PixelEngine::DimensionRange>(pyPixelEngine, "DimensionRange")
        .def(py::init<uint32_t, uint32_t, uint32_t>(),
            py::arg("begin"), py::arg("increment"),  py::arg("end"))
        .def_readonly("begin", &PixelEngine::DimensionRange::begin)
        .def_readonly("increment", &PixelEngine::DimensionRange::increment)
        .def_readonly("end", &PixelEngine::DimensionRange::end);
    
    py::class_<PixelEngine::BlockHeaderTemplate>(pyPixelEngine, "BlockHeaderTemplate")
        .def_readonly("bits_allocated", &PixelEngine::BlockHeaderTemplate::bitsAllocated)
        .def_readonly("bits_stored", &PixelEngine::BlockHeaderTemplate::bitsStored)
        .def_readonly("high_bit", &PixelEngine::BlockHeaderTemplate::highBit)
        .def_readonly("samples_per_pixel", &PixelEngine::BlockHeaderTemplate::samplesPerPixel)
        .def_readonly("compressor_id", &PixelEngine::BlockHeaderTemplate::compressorId)
        .def_readonly("dimensions_in_block", &PixelEngine::BlockHeaderTemplate::dimensionsInBlock)
        .def_readonly("dimension_ranges", &PixelEngine::BlockHeaderTemplate::dimensionRanges);

    py::class_<PixelEngine::BlockHeader>(pyPixelEngine, "BlockHeader")
        .def_readonly("coordinate", &PixelEngine::BlockHeader::coordinate)
        .def_readonly("template_id", &PixelEngine::BlockHeader::templateId)
        .def_readonly("offset", &PixelEngine::BlockHeader::offset)
        .def_readonly("size", &PixelEngine::BlockHeader::size);

    py::class_<PixelEngine::ClusterBlockHeader>(pyPixelEngine, "ClusterBlockHeader")
        .def_readonly("coordinate", &PixelEngine::ClusterBlockHeader::coordinate)
        .def_readonly("template_id", &PixelEngine::ClusterBlockHeader::templateId);

    py::class_<PixelEngine::ClusterHeaderTemplate>(pyPixelEngine, "ClusterHeaderTemplate")
        .def_readonly("dimensions_in_cluster", &PixelEngine::ClusterHeaderTemplate::dimensionsInCluster)
        .def_readonly("block_headers", &PixelEngine::ClusterHeaderTemplate::blockHeaders)
        .def_readonly("dimension_ranges", &PixelEngine::ClusterHeaderTemplate::dimensionRanges);

    py::class_<PixelEngine::ClusterHeader>(pyPixelEngine, "ClusterHeader")
        .def_readonly("coordinate", &PixelEngine::ClusterHeader::coordinate)
        .def_readonly("block_offsets", &PixelEngine::ClusterHeader::blockOffsets)
        .def_readonly("block_sizes", &PixelEngine::ClusterHeader::blockSizes)
        .def_readonly("template_id", &PixelEngine::ClusterHeader::templateId)
        .def_readonly("offset", &PixelEngine::ClusterHeader::offset)
        .def_readonly("size", &PixelEngine::ClusterHeader::size);

    py::class_<PixelEngine::SubImageHeader>(pyPixelEngine, "SubImageHeader")
        .def_readonly("cluster_header_templates", &PixelEngine::SubImageHeader::clusterHeaderTemplates)
        .def_readonly("cluster_sequence", &PixelEngine::SubImageHeader::clusterSequence)
        .def_readonly("block_header_templates", &PixelEngine::SubImageHeader::blockHeaderTemplates)
        .def_readonly("block_headers", &PixelEngine::SubImageHeader::blockHeaders)
        .def_readonly("cluster_headers", &PixelEngine::SubImageHeader::clusterHeaders);
        
    py::class_<PixelEngine::SubImage>(pyPixelEngine, "SubImage")
        .def("include_input_region", &PixelEngine::SubImage::includeInputRegion, py::arg("pos"), py::arg("size"), py::arg("level"))
        .def("preallocate_pixels", &PixelEngine::SubImage::preallocatePixels, py::arg("pos"), py::arg("size"), py::arg("level"))
        .def("put_block", [](PixelEngine::SubImage &self, py::buffer &buffer)
             { 
                py::buffer_info info = buffer.request();
                return self.putBlock(info.ptr, info.itemsize * info.size); },
            py::arg("buffer"))
        .def_property_readonly("has_display_view", &PixelEngine::SubImage::hasDisplayView)
        .def_property_readonly("source_view",
            py::overload_cast<>(&PixelEngine::SubImage::sourceView),
            py::return_value_policy::reference)
        .def_property_readonly("display_view",
            py::overload_cast<>(&PixelEngine::SubImage::displayView),
            py::return_value_policy::reference)
        .def("add_view", &PixelEngine::SubImage::addView, py::return_value_policy::reference)
        .def("block_size", &PixelEngine::SubImage::blockSize, py::arg("template_id") = 0)
        .def("block_pos", &PixelEngine::SubImage::blockPos, py::arg("block_ind"))
        .def("read_block", [](PixelEngine::SubImage &self, py::buffer &buffer)
             { 
                py::buffer_info info = buffer.request();
                return self.readBlock(info.ptr, info.itemsize * info.size); },
            py::arg("buffer"))
        .def("ordered_block_coordinates", &PixelEngine::SubImage::orderedBlockCoordinates)
        .def_property_readonly("image_type", &PixelEngine::SubImage::imageType)
        .def_property_readonly("pixel_transform", &PixelEngine::SubImage::pixelTransform)
        .def_property_readonly("quality_preset", &PixelEngine::SubImage::qualityPreset)
        .def_property_readonly("quality", &PixelEngine::SubImage::quality)
        .def_property_readonly("compressor", &PixelEngine::SubImage::compressor)
        .def_property_readonly("colorspace_transform", &PixelEngine::SubImage::colorspaceTransform)
        .def_property_readonly("num_tiles", &PixelEngine::SubImage::numTiles)
        .def_property("icc_profile",
            py::overload_cast<>(&PixelEngine::SubImage::iccProfile, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::SubImage::iccProfile))
        .def_property_readonly("icc_matrix", &PixelEngine::SubImage::iccMatrix)
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
        .def_property("lossy_image_compression",
            py::overload_cast<>(&PixelEngine::SubImage::lossyImageCompression, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::SubImage::lossyImageCompression))
        .def_property("lossy_image_compression_ratio",
            py::overload_cast<>(&PixelEngine::SubImage::lossyImageCompressionRatio, py::const_),
            py::overload_cast<double>(&PixelEngine::SubImage::lossyImageCompressionRatio))
        .def_property("lossy_image_compression_method",
            py::overload_cast<>(&PixelEngine::SubImage::lossyImageCompressionMethod, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::SubImage::lossyImageCompressionMethod))
        .def_property("color_linearity",
            py::overload_cast<>(&PixelEngine::SubImage::colorLinearity, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::SubImage::colorLinearity))
        .def("header", &PixelEngine::SubImage::header)
        .def("block_coordinate", &PixelEngine::SubImage::blockCoordinate, py::arg("block_index"))
        .def("block_index", &PixelEngine::SubImage::blockIndex, py::arg("block_coordinate"));

    py::class_<PixelEngine::ISyntaxFacade>(pyPixelEngine, "ISyntaxFacade")
        .def(
            "open", [](PixelEngine::ISyntaxFacade &self, std::filesystem::path const &url, std::string const &container_name, std::string const &mode, std::string const &cache_name)
            {
        std::ios_base::openmode open_mode = parse_isyntax_open_mode(mode);
        return self.open(url, container_name, open_mode, cache_name); },
            py::arg("url"),
            py::arg("container_name") = "",
            py::arg("mode") = "r",
            py::arg("cache_name") = "")
        .def(
            "open", [](PixelEngine::ISyntaxFacade &self, std::istream& stream, std::string const &container_name, std::string const &mode, std::string const &cache_name)
            {
        std::ios_base::openmode open_mode = parse_isyntax_open_mode(mode);
        return self.open(static_cast<std::iostream*>(&stream), container_name, open_mode, cache_name); },
            py::arg("stream"),
            py::arg("container_name") = "",
            py::arg("mode") = "r",
            py::arg("cache_name") = "")
        .def("add_sub_image", &PixelEngine::ISyntaxFacade::addSubImage)
        .def_property_readonly("num_images", &PixelEngine::ISyntaxFacade::numImages)
        .def("__getitem__",
            py::overload_cast<size_t>(&PixelEngine::ISyntaxFacade::operator[]),
            py::arg("index"),
            py::return_value_policy::reference)
        .def("__getitem__",
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::operator[]),
            py::arg("type"),
            py::return_value_policy::reference)
        .def_property_readonly("isyntax_file_version", &PixelEngine::ISyntaxFacade::iSyntaxFileVersion)
        .def_property_readonly("id", &PixelEngine::ISyntaxFacade::id)
        .def_property("barcode",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::barcode, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::barcode))
        .def_property("scanner_calibration_status",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::scannerCalibrationStatus, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::scannerCalibrationStatus))
        .def_property("software_versions",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::softwareVersions, py::const_),
            py::overload_cast<std::vector<std::string> const&>(&PixelEngine::ISyntaxFacade::softwareVersions))
        .def_property("derivation_description",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::derivationDescription, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::derivationDescription))
        .def_property("acquisition_datetime",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::acquisitionDateTime, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::acquisitionDateTime))
        .def_property("manufacturer",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::manufacturer, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::manufacturer))
        .def_property("model_name",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::modelName, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::modelName))
        .def_property("device_serial_number",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::deviceSerialNumber, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::deviceSerialNumber))
        .def_property("scanner_rack_number",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::scannerRackNumber, py::const_),
            py::overload_cast<uint16_t>(&PixelEngine::ISyntaxFacade::scannerRackNumber))
        .def_property("scanner_slot_number",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::scannerSlotNumber, py::const_),
            py::overload_cast<uint16_t>(&PixelEngine::ISyntaxFacade::scannerSlotNumber))
        .def_property("scanner_operator_id",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::scannerOperatorId, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::ISyntaxFacade::scannerOperatorId))
        .def_property("scanner_rack_priority",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::scannerRackPriority, py::const_),
            py::overload_cast<uint16_t>(&PixelEngine::ISyntaxFacade::scannerRackPriority))
        .def_property("date_of_last_calibration",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::dateOfLastCalibration, py::const_),
            py::overload_cast<std::vector<std::string> const&>(&PixelEngine::ISyntaxFacade::dateOfLastCalibration))
        .def_property("time_of_last_calibration",
            py::overload_cast<>(&PixelEngine::ISyntaxFacade::timeOfLastCalibration, py::const_),
            py::overload_cast<std::vector<std::string> const&>(&PixelEngine::ISyntaxFacade::timeOfLastCalibration))
        .def_property_readonly("is_philips", &PixelEngine::ISyntaxFacade::isPhilips)
        .def_property_readonly("is_hamamatsu", &PixelEngine::ISyntaxFacade::isHamamatsu)
        .def_property_readonly("is_UFS", &PixelEngine::ISyntaxFacade::isUFS)
        .def_property_readonly("is_UFSb", &PixelEngine::ISyntaxFacade::isUFSb)
        .def_property_readonly("is_UVS", &PixelEngine::ISyntaxFacade::isUVS)
        .def("close", &PixelEngine::ISyntaxFacade::close);

    pyPixelEngine.def(py::init<>())
        .def(py::init<RenderBackend &, RenderContext &>(),
            py::arg("render_backend"), py::arg("render_context"))
        .def_property_readonly_static("version", [](py::object)
                                      { return PixelEngine::version(); })
        .def("__getitem__", &PixelEngine::operator[], py::arg("name"),
            py::return_value_policy::reference)
        .def("containers", &PixelEngine::containers)
        .def("container_version", &PixelEngine::containerVersion, py::arg("container"))
        .def("compressors", &PixelEngine::compressors)
        .def("pixel_transforms", &PixelEngine::pixelTransforms)
        .def("block_sizes", &PixelEngine::blockSizes, py::arg("pixel_transform"))
        .def("colorspace_transforms", &PixelEngine::colorspaceTransforms)
        .def("quality_presets", &PixelEngine::qualityPresets)
        .def("supported_filters", &PixelEngine::supportedFilters)
        .def("wait_all", &PixelEngine::waitAll, py::arg("regions"))
        .def("wait_any",
            py::overload_cast<>(&PixelEngine::waitAny),
            py::return_value_policy::reference)
        .def("wait_any", 
            py::overload_cast<std::list<std::shared_ptr<PixelEngine::Region>> const&>(&PixelEngine::waitAny),
            py::arg("regions"),
            py::return_value_policy::reference)
        .def("clear_render_target", &PixelEngine::clearRenderTarget,
            py::arg("color"),
            py::arg("target") = 0)
        .def("clear_render_cache", &PixelEngine::clearRenderCache)
        .def("clear_render_buffers", &PixelEngine::clearRenderBuffers)
        .def_property("network_timeout",
            py::overload_cast<>(&PixelEngine::networkTimeout, py::const_),
            py::overload_cast<size_t>(&PixelEngine::networkTimeout))
        .def_property("client_certificates",
            py::overload_cast<>(&PixelEngine::clientCertificates, py::const_),
            [](PixelEngine &self, std::tuple<std::string, std::string, std::string> const &certs)
            {
                auto const& [cert, key, password] = certs;
                return self.clientCertificates(cert, key, password); 
            })
        .def_property("certificates",
            py::overload_cast<>(&PixelEngine::certificates, py::const_),
            py::overload_cast<std::string const&>(&PixelEngine::certificates))
        .def_property("rate_limiting",
            py::overload_cast<>(&PixelEngine::rateLimiting, py::const_),
            py::overload_cast<bool>(&PixelEngine::rateLimiting));
}