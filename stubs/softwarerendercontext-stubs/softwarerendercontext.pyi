import pixelengine

from typing import overload

class SoftwareRenderContext(pixelengine.RenderContext):
    number_of_worker_threads: int
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, width: int, height: int) -> None: ...
