from setuptools import Extension, setup
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "placeback_subvolume_cyt",
        ["placeback_subvolume_cyt.pyx"],
        libraries=["m"],
        extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
        extra_link_args=['-fopenmp']
    )
]

setup(
    name='placeback_subvolume_cyt',
    ext_modules=cythonize(ext_modules),
)