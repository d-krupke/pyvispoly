# pyvispoly

This package provides a Python interface to the [CGAL](https://www.cgal.org/) library for computing visibility polygons.
It is exact, but not necessarily efficient.

**Motivation:** This package was developed in the context of implementing an
exact solver for the dispersive Art Gallery Problem.
Because the problem is difficult and only reasonably small instances can be solved,
the focus of this package is on correctness and not on efficiency.

## Installation

You can install the package using pip:

```bash
pip install --verbose .
```

During installation, it will download and install CGAL and its dependencies. This may take a while.

## Challenges of visibility polygon computation

Computing and working with visibility polygons requires exact arithmetic, and will
otherwise frequently lead to incorrect results. Luckily, the CGAL library provides
a robust implementation of visibility polygon computation, which is used by this package.
However, you need to convert your coordinates to CGAL's exact number type.
If you forget to use the correct types, you may get exceptions.
This library does not support duck typing, and will not convert your coordinates
to the correct type. While automatic conversion would be possible, it could hide
errors in your code. When working with geometry, you have to be exact and not work with
Pythons inexact floating point numbers.

## Usage

```python
# import elements from pyvispoly
from pyvispoly import FieldNumber, Point, Polygon, PolygonWithHoles, VisibilityPolygonCalculator, plot_polygon
import matplotlib.pyplot as plt

# Polygon with holes
poly = PolygonWithHoles([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)],
                         [[Point(0.25, 0.25), Point(0.75, 0.25), Point(0.75, 0.75), Point(0.25, 0.75)][::-1]])

# compute visibility polygon
visp_poly_calc = VisibilityPolygonCalculator(poly)
vis_poly = visp_poly_calc.compute_visibility_polygon(Point(0.2, 0.0))

fig, ax = plt.subplots()
ax.set_aspect('equal')
plt.title("Visibility Polygon")
plot_polygon(poly, ax=ax, color="lightgrey")
plot_polygon(vis_poly, ax=ax, color='red', alpha=0.5)
plt.plot([0.2], [0.0], 'x', color='black')
plt.show()
```

![Visibility Polygon](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/visibility_polygon.png?raw=true)

## License

This library statically links to parts of CGAL, which are licensed under GPL, making this library also GPL.
If you have a special (commercial) license for CGAL, this library (without the CGAL component)
is under the more generous MIT license (as most of my code).
