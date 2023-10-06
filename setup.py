from skbuild_conan import setup
from setuptools import find_packages

setup(  # https://scikit-build.readthedocs.io/en/latest/usage.html#setup-options
    name="pyvispoly",
    version="0.1.1",
    author="Dominik Krupke",
    license="LICENSE",
    description="README.md",
    long_description="README.md",
    packages=find_packages("src"),  # Include all packages in `./src`.
    package_dir={"": "src"},  # The root for our python package is in `./src`.
    python_requires=">=3.7",  # lowest python version supported.
    install_requires=[
        "matplotlib>=3.4.2",
    ],  # Python Dependencies
    conan_requirements=["fmt/[>=10.0.0]", "cgal/[>=5.6]"],  # C++ Dependencies
    cmake_minimum_required_version="3.23",
)
