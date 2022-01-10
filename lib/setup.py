from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "placeback_subvolume_cyt",
        ["placeback_subvolume_cyt.pyx"],
        libraries=["m"],
        include_dirs=[numpy.get_include()],
        extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
        extra_link_args=['-fopenmp']
    )
]

setup(
    name='placeback_subvolume_cyt',
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(ext_modules),
)
