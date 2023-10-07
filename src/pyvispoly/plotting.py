"""
Some plotting utilities for polygons.

2023, Dominik Krupke, TU Braunschweig
"""
import typing

import matplotlib.pyplot as plt

from ._cgal_bindings import Polygon, PolygonWithHoles


def plot_polygon(
    polygon: typing.Union[Polygon, PolygonWithHoles], color="blue", ax=None, **kwargs
):
    if ax is None:
        ax = plt.gca()
    if type(polygon) == PolygonWithHoles:
        plot_polygon(polygon.outer_boundary(), color=color, ax=ax, **kwargs)
        for hole in polygon.holes():
            plot_polygon(hole, color="white")
    else:
        x = [float(p.x()) for p in polygon.boundary()]
        y = [float(p.y()) for p in polygon.boundary()]
        ax.fill(x, y, color=color, **kwargs)
