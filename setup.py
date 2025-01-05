from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "2.0"

ext_modules = [
    Pybind11Extension("pixelengine",
        ["src/pixelengine.cpp"],
        language='c++',
        cxx_std=14,
        libraries=['pixelengine']
        ),
    Pybind11Extension("softwarerenderbackend",
        ["src/softwarerenderbackend.cpp"],
        language='c++',
        cxx_std=14,
        libraries=['softwarerenderbackend']
        ),
    Pybind11Extension("softwarerendercontext",
        ["src/softwarerendercontext.cpp"],
        language='c++',
        cxx_std=14,
        libraries=['softwarerendercontext']
        ),
]

setup(
    name="philips-pathologysdk-python",
    version=__version__,
    author="Stefano Malacrino",
    url="https://github.com/stefano-malacrino/philips-pathologysdk-python",
    description="Python bindings for the Philips Pathology SDK",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.6",
)
