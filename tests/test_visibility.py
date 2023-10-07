from pyvispoly import (
    FieldNumber as FT,
)
from pyvispoly import (
    Point,
    Polygon,
    PolygonWithHoles,
    VisibilityPolygonCalculator,
)


def test_contains():
    boundary = [
        Point(FT(0), FT(0)),
        Point(FT(0), FT(1)),
        Point(FT(1), FT(1)),
        Point(FT(1), FT(0)),
    ]
    polygon = Polygon(boundary[::-1])
    assert polygon.contains(Point(FT(0.5), FT(0.5)))
    assert not polygon.contains(Point(FT(1.5), FT(0.5)))
    print(str(polygon.area()))
    assert polygon.area() > FT(0)
    polygon_with_holes = PolygonWithHoles(polygon, [])
    vispoly = VisibilityPolygonCalculator(polygon_with_holes)
    assert not vispoly.is_feasible_query_point(Point(FT(1.5), FT(0.5)))
    assert vispoly.is_feasible_query_point(Point(FT(0.5), FT(0.5)))
    vis = vispoly.compute_visibility_polygon(Point(FT(0.5), FT(0.5)))
    assert vis.area() == FT(1)
