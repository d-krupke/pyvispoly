from pyvispoly import (
    FieldNumber,
    Point,
    Polygon,
    PolygonWithHoles,
    VisibilityPolygonCalculator,
)


def test_contains():
    boundary = [
        Point(FieldNumber(0), FieldNumber(0)),
        Point(FieldNumber(0), FieldNumber(1)),
        Point(FieldNumber(1), FieldNumber(1)),
        Point(FieldNumber(1), FieldNumber(0)),
    ]
    polygon = Polygon(boundary[::-1])
    assert polygon.contains(Point(FieldNumber(0.5), FieldNumber(0.5))) == True
    assert polygon.contains(Point(FieldNumber(1.5), FieldNumber(0.5))) == False
    print(str(polygon.area()))
    assert polygon.area() > FieldNumber(0)
    polygon_with_holes = PolygonWithHoles(polygon, [])
    vispoly = VisibilityPolygonCalculator(polygon_with_holes)
    assert (
        vispoly.is_feasible_query_point(Point(FieldNumber(1.5), FieldNumber(0.5)))
        == False
    )
    assert (
        vispoly.is_feasible_query_point(Point(FieldNumber(0.5), FieldNumber(0.5)))
        == True
    )
    vis = vispoly.compute_visibility_polygon(Point(FieldNumber(0.5), FieldNumber(0.5)))
    assert vis.area() == FieldNumber(1)
