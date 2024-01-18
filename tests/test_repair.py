from pyvispoly import (
    FieldNumber as FT,
)
from pyvispoly import (
    Point,
    Polygon,
    PolygonWithHoles,
    repair,
)


def test_repair():
    boundary = [
        Point(FT(0), FT(0)),
        Point(FT(0), FT(1)),
        Point(FT(1), FT(1)),
        Point(FT(1), FT(0)),
    ]
    polygon = Polygon(boundary[::-1])
    poly_with_holes = PolygonWithHoles(polygon, [])
    poly_ = repair(poly_with_holes)[0]
    assert poly_.area() == 1.0


def test_repair_with_hole():
    boundary = [
        Point(FT(0), FT(0)),
        Point(FT(0), FT(1)),
        Point(FT(1), FT(1)),
        Point(FT(1), FT(0)),
    ]
    hole = [
        Point(FT(0.25), FT(0.25)),
        Point(FT(0.25), FT(0.75)),
        Point(FT(0.75), FT(0.75)),
        Point(FT(0.75), FT(0.25)),
    ]
    polygon = Polygon(boundary[::-1])
    poly_with_holes = PolygonWithHoles(polygon, [Polygon(hole)])
    poly_ = repair(poly_with_holes)[0]
    assert all(hole.area() >= FT(0) for hole in poly_.holes())
    assert poly_.area() == 0.75
