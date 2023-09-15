# pyvispoly

This package provides a Python interface to the [CGAL](https://www.cgal.org/) library for computing visibility polygons.
It is exact, but not necessarily efficient.

> :warning: While the project code itself is free open source software (MIT), it depends on CGAL which may require a license for commercial use. See [CGAL license](https://www.cgal.org/license.html) for more information.

## Installation

You can install the package using pip:

```bash
pip install --verbose .
```

During installation, it will download and install CGAL and its dependencies. This may take a while.

## Usage

```python
from pyvispoly import VisibilityPolygon, Point, PolygonWithHoles
from pyvispoly.plotting import plot_polygon

polygon = PolygonWithHoles(
    outer_boundary=[Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)],
    holes=[
        [Point(0.25, 0.25), Point(0.75, 0.25), Point(0.75, 0.75), Point(0.25, 0.75)]
    ],
)

vis_calculator = VisibilityPolygon(polygon)
vis_polygon = vis_calculator.compute_visibility_polygon(Point(0.5, 0.5))

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
plot_polygon(polygon, ax=ax, alpha=0.2, color="blue")
plot_polygon(vis_polygon, ax=ax, alpha=0.5, color="red")
ax.plot(0.5, 0.5, "o", color="black")
ax.set_aspect('equal', 'box')
plt.show()
```
