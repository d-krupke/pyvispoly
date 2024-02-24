# pyvispoly

This package provides a Python interface to the [CGAL](https://www.cgal.org/)
library for computing visibility polygons. It is exact, but not necessarily
efficient.

**Motivation:** This package was developed in the context of implementing an
exact solver for the dispersive Art Gallery Problem
([see here](https://github.com/d-krupke/dispersive_agp_solver)). Because the
problem is difficult and only reasonably small instances can be solved, the
focus of this package is on correctness and not on efficiency.

![Visibility Polygon](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/visibility_polygon.png?raw=true)

## Installation

You can install the package using pip:

```bash
pip install --verbose .
```

or directly via PyPI:

```bash
pip install --verbose pyvispoly
```

During installation, it will download and install CGAL and its dependencies.
This may take a while. This requires a C++-compiler to be installed on your
system.

> :info: On installation problems, please first check out
> [skbuild-conan's troubleshooting guide](https://github.com/d-krupke/skbuild-conan#common-problems).
> It is likely that the automatic C++-compilation fails.

## Addressing the Complexity of Visibility Polygon Computation

Computing visibility polygons demands precise arithmetic; imprecise calculations
often yield incorrect results. This package leverages the robustness of the CGAL
library for accurate visibility polygon computation. Note that it is essential
to convert your coordinates to CGAL's exact number type. Failure to use the
correct types may result in exceptions.

This package explicitly avoids duck typing and automatic type conversion. While
these features could simplify usage, they risk obscuring underlying errors in
your code. In geometric computations, precision is paramount; thus, Python's
native floating-point numbers are unsuitable.

## Usage

```python
# import elements from pyvispoly
from pyvispoly import (
    FieldNumber,
    Point,
    Polygon,
    PolygonWithHoles,
    VisibilityPolygonCalculator,
    plot_polygon,
)
import matplotlib.pyplot as plt
```

```python
# exact number type
a = FieldNumber(100)
b = FieldNumber(200)
c = a + b
assert isinstance(c, FieldNumber)
print("c:", c)
print("float(c):", float(c))
assert float(c) == 300.0
```

    c: 300.000000
    float(c): 300.0

```python
# point type
p1 = Point(FieldNumber(0), FieldNumber(0))
assert p1.x() == FieldNumber(0) and p1.y() == FieldNumber(0)
p2 = Point(2, 3)
p3 = Point(4.0, 5.0)
print("p1:", p1)
print("p2:", p2)
print("p3:", p3)
```

    p1: (0, 0)
    p2: (2, 3)
    p3: (4, 5)

```python
# Polygon
poly1 = Polygon([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)])
assert poly1.is_simple()
assert poly1.area() == FieldNumber(1)
assert float(poly1.area()) == 1.0
```

```python
# The ccw order of the points is important
poly1 = Polygon([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)][::-1])
assert poly1.is_simple()
assert poly1.area() == FieldNumber(-1)
assert float(poly1.area()) == -1.0
```

```python
# Polygon with holes
poly2 = PolygonWithHoles(
    [Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)],
    [
        [Point(0.25, 0.25), Point(0.75, 0.25), Point(0.75, 0.75), Point(0.25, 0.75)][
            ::-1
        ]
    ],
)
# boundary is counter-clockwise
assert poly2.outer_boundary().area() >= FieldNumber(0)
for hole in poly2.holes():
    # holes are clockwise
    assert hole.area() <= FieldNumber(0)

fig, ax = plt.subplots()
ax.set_aspect("equal")
plot_polygon(poly2, ax=ax)
plt.show()
```

![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_5_0.png?raw=true)

```python
# set operations on polygons
poly3 = PolygonWithHoles([Point(0, 0), Point(0.1, 0), Point(0.1, 0.1), Point(0, 0.1)])
poly4_list = poly2.difference(poly3)
assert len(poly4_list) == 1, "Expect a single polygon"
poly4 = poly4_list[0]
fig, ax = plt.subplots()
ax.set_aspect("equal")
plot_polygon(poly4, ax=ax)
plt.show()
```

![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_6_0.png?raw=true)

```python
# compute visibility polygon
visp_poly_calc = VisibilityPolygonCalculator(poly4)
vis_poly = visp_poly_calc.compute_visibility_polygon(Point(0.2, 0.0))

fig, ax = plt.subplots()
ax.set_aspect("equal")
plt.title("Visibility Polygon")
plot_polygon(poly4, ax=ax, color="lightgrey")
plot_polygon(vis_poly, ax=ax, color="red", alpha=0.5)
plt.plot([0.2], [0.0], "x", color="black")
plt.savefig("../docs/figures/visibility_polygon.png")
plt.show()
```

![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_7_0.png?raw=true)

```python
# getting some interior points
fig, ax = plt.subplots()
ax.set_aspect("equal")
plt.title("Sample points in the interior")
plot_polygon(poly4, ax=ax, color="lightgrey")
plot_polygon(vis_poly, ax=ax, color="red", alpha=0.5)
plt.plot([0.2], [0.0], "x", color="black")
interior_points = vis_poly.interior_sample_points()
for p in interior_points:
    assert vis_poly.contains(p), "all points should be inside the visibility polygon"
plt.plot(
    [p.x() for p in interior_points],
    [p.y() for p in interior_points],
    "+",
    color="black",
)
plt.show()
```

![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_8_0.png?raw=true)

> :warning: The library may segfault if you pass bad polygons. For example, if
> holes intersect the outer boundary, the library may crash. One could probably
> catch these errors, but I have not implemented this yet (feel free to do so
> and submit a pull request).

## License

This library incorporates components from the CGAL library through static
linking, which are subject to the terms of the GNU General Public License (GPL).
Consequently, the use and distribution of this library are also governed by the
GPL.

However, if you possess a specialized (commercial) license for CGAL (or replace
CGAL), then the portions of this library excluding the CGAL component are
available under the more permissive MIT license, consistent with the licensing
of my other code.

Please ensure you understand the implications of these licenses before using or
distributing this library.

## Changelog

- **0.3.0** Added `do_intersect`.
- **0.2.0** Made the library more robust and added a `repair`-method as the
  boolean operations are not always guaranteed to produce valid polygons. Note
  that this is just kind of a hack and may not always work.
- **0.1.3:** (2023-11-29) Workaround to bug in CGAL.
