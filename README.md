# pyvispoly

This package provides a Python interface to the [CGAL](https://www.cgal.org/) library for computing visibility polygons.
It is exact, but not necessarily efficient.

**Motivation:** This package was developed in the context of implementing an
exact solver for the dispersive Art Gallery Problem ([see here](https://github.com/d-krupke/dispersive_agp_solver)).
Because the problem is difficult and only reasonably small instances can be solved,
the focus of this package is on correctness and not on efficiency.

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
```


```python
# exact number type
a = FieldNumber(100)
b = FieldNumber(200)
c = a+b
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
poly2 = PolygonWithHoles([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)],
                         [[Point(0.25, 0.25), Point(0.75, 0.25), Point(0.75, 0.75), Point(0.25, 0.75)][::-1]])
# boundary is counter-clockwise
assert poly2.outer_boundary().area() >= FieldNumber(0)
for hole in poly2.holes():
    # holes are clockwise
    assert hole.area() <= FieldNumber(0)

fig, ax = plt.subplots()
ax.set_aspect('equal')
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
ax.set_aspect('equal')
plot_polygon(poly4, ax=ax)
plt.show()
```


    
![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_6_0.png?raw=true)


```python
# compute visibility polygon
visp_poly_calc = VisibilityPolygonCalculator(poly4)
vis_poly = visp_poly_calc.compute_visibility_polygon(Point(0.2, 0.0))

fig, ax = plt.subplots()
ax.set_aspect('equal')
plt.title("Visibility Polygon")
plot_polygon(poly4, ax=ax, color="lightgrey")
plot_polygon(vis_poly, ax=ax, color='red', alpha=0.5)
plt.plot([0.2], [0.0], 'x', color='black')
plt.savefig("../docs/figures/visibility_polygon.png")
plt.show()
```


![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_7_0.png?raw=true)



```python
# getting some interior points
fig, ax = plt.subplots()
ax.set_aspect('equal')
plt.title("Sample points in the interior")
plot_polygon(poly4, ax=ax, color="lightgrey")
plot_polygon(vis_poly, ax=ax, color='red', alpha=0.5)
plt.plot([0.2], [0.0], 'x', color='black')
interior_points = vis_poly.interior_sample_points()
for p in interior_points:
    assert vis_poly.contains(p), "all points should be inside the visibility polygon"
plt.plot([p.x() for p in interior_points], [p.y() for p in interior_points], '+', color='black')
plt.show()
```

![png](https://github.com/d-krupke/pyvispoly/blob/main/docs/figures/simple_example_8_0.png?raw=true)


## License

This library statically links to parts of CGAL, which are licensed under GPL, making this library also GPL.
If you have a special (commercial) license for CGAL, this library (without the CGAL component)
is under the more generous MIT license (as most of my code).
