import os
import setuptools
import sys

from pybind11.setup_helpers import Pybind11Extension, build_ext

setuptools.setup(
    ext_modules=[
        Pybind11Extension(
            "fastpyfqmr",
            ["pyfqmr/fastpyfqmr.cpp"],
            extra_compile_args = ['-std=c++20'],
            language="c++",
        ),
    ],
)



