from setuptools import find_packages
from skbuild_conan import setup


def readme():
    """
    :return: Content of README.md
    """
    with open("README.md") as file:
        return file.read()


setup(
    # https://scikit-build.readthedocs.io/en/latest/usage.html#setup-options
    # https://github.com/d-krupke/skbuild-conan#usage
    name="pyvispoly",
    version="0.3.0",
    author="Dominik Krupke",
    license="LICENSE",
    description="CGAL's visibility polygons in Python",
    long_description=readme(),
    long_description_content_type="text/markdown",
    packages=find_packages("src"),  # Include all packages in `./src`.
    package_dir={"": "src"},  # The root for our python package is in `./src`.
    python_requires=">=3.7",  # lowest python version supported.
    install_requires=[
        "matplotlib>=3.4.2",
    ],  # Python Dependencies
    conan_requirements=["fmt/[>=10.0.0]", "cgal/[>=5.6]"],  # C++ Dependencies
    cmake_minimum_required_version="3.23",
)
