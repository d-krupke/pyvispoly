
// pybind11
#include <pybind11/operators.h> // to define operator overloading
#include <pybind11/pybind11.h>  // basic pybind11 functionality
#include <pybind11/stl.h>       // automatic conversion of vectors
// cgal
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
//#include <CGAL/Arr_naive_point_location.h>
// fmt
#include <fmt/core.h>

// Define CGAL types
using Kernel = CGAL::Epeck; // Exact Predicates Exact Constructions Kernel
using Point = CGAL::Point_2<Kernel>;
using Polygon2WithHoles = CGAL::Polygon_with_holes_2<Kernel>;
using Polygon2 = CGAL::Polygon_2<Kernel>;

using Segment2 = Kernel::Segment_2;
using Traits_2 = CGAL::Arr_segment_traits_2<Kernel>;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;
using Halfedge_const_handle = Arrangement_2::Halfedge_const_handle;
using Face_handle = Arrangement_2::Face_handle;

// Define the used visibility class
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2> TEV;

class VisibilityPolygonCalculator
{
public:
  VisibilityPolygonCalculator(Polygon2WithHoles &poly)
  {
    if (!poly.outer_boundary().is_simple())
    {
      throw std::runtime_error("Polygon is not simple");
    }
    polygon = poly;
    std::vector<Segment2> segments;
    if (poly.outer_boundary().area() < 0)
    {
      throw std::runtime_error("Polygon is not counterclockwise oriented");
    }
    for (const auto e : poly.outer_boundary().edges())
    {
      auto s = e.source();
      auto t = e.target();
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
    for (const auto &hole : poly.holes())
    {
      if (hole.area() > 0)
      {
        throw std::runtime_error("Hole is not clockwise oriented");
      }
      for (const auto e : hole.edges())
      {
        auto s = e.source();
        auto t = e.target();
        auto seg = Segment2(s, t);
        segments.push_back(seg);
      }
    }
    CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());
    auto face = env.unbounded_face();
    if (face->number_of_holes() != 1 || !face->is_unbounded())
    {
      throw std::runtime_error("Bad arrangement");
    }
    interior_face = face->holes_begin()->ptr()->twin()->face();
    if (interior_face->is_unbounded())
    {
      throw std::runtime_error("Bad arrangement");
    }
  }



  bool is_feasible_query_point(const Point &query_point)
  {
    if (polygon.outer_boundary().oriented_side(query_point) == CGAL::ON_NEGATIVE_SIDE)
    {
      return false;
    }
    for (const auto &hole : polygon.holes())
    {
      if (hole.oriented_side(query_point) == CGAL::ON_POSITIVE_SIDE)
      {
        return false;
      }
    }
    return true;
  }

  Polygon2 compute_visibility_polygon(const Point &query_point)
  {
    if (!is_feasible_query_point(query_point))
    {
      throw std::runtime_error("Query point not feasible");
    }
    Arrangement_2 output_arr;
    TEV tev(env);
    Face_handle fh = tev.compute_visibility(query_point, interior_face, output_arr);
    Polygon2 poly;
    std::vector<Point> points;
    auto e = fh->outer_ccb();
    points.push_back(e->source()->point());
    ++e;
    while (e != fh->outer_ccb())
    {
      points.push_back(e->source()->point());
      ++e;
    }
    return Polygon2(points.begin(), points.end());
  }

  Polygon2WithHoles polygon;
  Arrangement_2 env;
  Face_handle interior_face;
};

// Getting this name right is important! It has to equal the name in the
// CMakeLists.txt.
PYBIND11_MODULE(_cgal_bindings, m)
{
  namespace py = pybind11;
  m.doc() = "Example of PyBind11 and CGAL."; // optional module docstring

  // Exact numbers
  py::class_<Kernel::FT>(m, "FieldNumber",
                         "A container for exact numbers in CGAL.")
      .def(py::init<long>())
      .def(py::init<double>())
      .def(py::self / Kernel::FT())
      .def(py::self + Kernel::FT())
      .def(py::self * Kernel::FT())
      .def(py::self == Kernel::FT())
      .def(py::self < Kernel::FT())
      .def(py::self > Kernel::FT())
      .def(py::self <= Kernel::FT())
      .def(py::self >= Kernel::FT())
      .def("__float__", &CGAL::to_double<Kernel::FT>)
      .def("__str__", [](const Kernel::FT &x)
           { return std::to_string(CGAL::to_double(x)); });

  // Points
  py::class_<Point>(m, "Point", "A point in CGAL.")
      .def(py::init<long, long>())
      .def(py::init<double, double>())
      .def(py::init<Kernel::FT, Kernel::FT>())
      .def("x", [](const Point &p)
           { return p.x(); })
      .def("y", [](const Point &p)
           { return p.y(); })
      .def(py::self == Point())
      .def("__str__", [](const Point &p)
           { return fmt::format("({}, {})", CGAL::to_double(p.x()),
                                CGAL::to_double(p.y())); });

  // Polygons
  py::class_<Polygon2>(m, "Polygon", "A simple polygon in CGAL.")
      .def(py::init<>())
      .def(py::init([](const std::vector<Point> &vertices)
                    { return std::make_unique<Polygon2>(vertices.begin(), vertices.end()); }))
      .def("boundary",
           [](const Polygon2 &poly)
           {
             std::vector<Point> points;
             std::copy(poly.begin(), poly.end(), std::back_inserter(points));
             return points;
           })
      .def("is_simple", &Polygon2::is_simple)
      .def("contains", [](const Polygon2 &self, const Point &p)
           { return self.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE; })
      .def("area", [](const Polygon2 &poly)
           { return poly.area(); });

  py::class_<Polygon2WithHoles>(m, "PolygonWithHoles",
                                "A polygon with holes in CGAL.")
      .def(py::init(
          [](const Polygon2 &outer, const std::vector<Polygon2> &holes)
          {
            return new Polygon2WithHoles(outer, holes.begin(), holes.end());
          }))
      .def("outer_boundary",
           [](const Polygon2WithHoles &poly)
           { return poly.outer_boundary(); })
      .def("holes", [](const Polygon2WithHoles &poly)
           {
        std::vector<Polygon2> holes;
        std::copy(poly.holes_begin(), poly.holes_end(),
                  std::back_inserter(holes));
        return holes; });

  py::class_<VisibilityPolygonCalculator>(m, "VisibilityPolygonCalculator",
                                          "A class to compute visibility polygons.")
      .def(py::init<Polygon2WithHoles &>())
      .def("compute_visibility_polygon", &VisibilityPolygonCalculator::compute_visibility_polygon)
      .def("is_feasible_query_point", &VisibilityPolygonCalculator::is_feasible_query_point);
}
