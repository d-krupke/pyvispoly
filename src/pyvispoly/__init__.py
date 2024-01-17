"""
This package provides the visibility polygons of CGAL as a
python package.

2023, Dominik Krupke, TU Braunschweig
"""
# ._cgal_bindings will only exist after compilation.
from ._cgal_bindings import (
    FieldNumber,
    Point,
    Polygon,
    PolygonWithHoles,
    VisibilityPolygonCalculator,
    repair,
)
from .plotting import plot_polygon

__all__ = [
    "FieldNumber",
    "Point",
    "Polygon",
    "PolygonWithHoles",
    "VisibilityPolygonCalculator",
    "plot_polygon",
    "repair",
]
