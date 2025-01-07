import pixelengine

class SoftwareRenderBackend(pixelengine.RenderBackend):
    def __init__(
        self,
        image_format_type: pixelengine.RenderBackend.ImageFormatType = pixelengine.RenderBackend.ImageFormatType.RGBA,
        tile_cache_size_bytes: int = 256000000,
    ) -> None: ...
