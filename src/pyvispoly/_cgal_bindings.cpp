/**
 * @file _cgal_bindings.cpp
 * @author Dominik Krupke (d.krupke@tu-bs.de)
 * @brief This file contains the pybind11 bindings for the CGAL library, which
 *      is used to compute visibility polygons. We need CGAL with it exact
 *      computations to ensure that the visibility polygons are computed
 *     correctly.
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023
 *
 */

// pybind11
#include <pybind11/operators.h> // to define operator overloading
#include <pybind11/pybind11.h>  // basic pybind11 functionality
#include <pybind11/stl.h>       // automatic conversion of vectors
// cgal
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/create_straight_skeleton_2.h>
//  fmt
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
using PointLocation = CGAL::Arr_naive_point_location<Arrangement_2>;
// Define the used visibility class
typedef CGAL::Triangular_expansion_visibility_2<
    Arrangement_2, /*Regularization*/ CGAL::Tag_true>
    TEV;
// using TEV = CGAL::Rotational_sweep_visibility_2<Arrangement_2>;

std::vector<Point> get_sample_points(Polygon2 &poly) {
  if (!poly.is_simple()) {
    throw std::runtime_error("Polygon is not simple.");
  }
  if (poly.area() <= 0) {
    throw std::runtime_error("Polygon has non-positive area. Cannot compute "
                             "sample points.");
  }
  std::vector<Point> points;
  // compute skeleton
  auto skeleton = CGAL::create_interior_straight_skeleton_2(
      poly.vertices_begin(), poly.vertices_end(), Kernel());
  // use skeleton vertices as sample points
  for (auto v = skeleton->vertices_begin(); v != skeleton->vertices_end();
       ++v) {
    auto p = v->point();
    // check if point is on boundary
    if (poly.bounded_side(p) == CGAL::ON_BOUNDARY) {
      // we only want internal points
      continue;
    }
    points.push_back(v->point());
  }
  return points;
}

void _check_polygon(Polygon2WithHoles &poly) {
  // Check if polygon is valid
  if (!poly.outer_boundary().is_simple()) {
    throw std::runtime_error("Polygon is not simple.");
  }
  if (poly.outer_boundary().area() <= 0) {
    throw std::runtime_error("Polygon has non-positive area. Cannot compute "
                             "sample points.");
  }
  if (std::any_of(poly.holes_begin(), poly.holes_end(),
                  [](const Polygon2 &hole) { return !hole.is_simple(); })) {
    throw std::runtime_error("Hole is not simple.");
  }
  if (std::any_of(poly.holes_begin(), poly.holes_end(),
                  [](const Polygon2 &hole) { return hole.area() >= 0; })) {
    throw std::runtime_error("Hole has non-negative area.");
  }
}

std::vector<Point> get_sample_points_with_holes(Polygon2WithHoles &poly) {
  // if no holes, use the other function
  if (poly.number_of_holes() == 0) {
    return get_sample_points(poly.outer_boundary());
  }
  _check_polygon(poly);
  std::vector<Point> points;
  // compute skeleton
  auto skeleton = CGAL::create_interior_straight_skeleton_2(
      poly.outer_boundary().vertices_begin(),
      poly.outer_boundary().vertices_end(), poly.holes_begin(),
      poly.holes_end(), Kernel());
  // use skeleton vertices as sample points
  for (auto v = skeleton->vertices_begin(); v != skeleton->vertices_end();
       ++v) {
    auto p = v->point();

    // check if point is on boundary
    if (poly.outer_boundary().bounded_side(p) == CGAL::ON_BOUNDARY) {
      // we only want internal points
      continue;
    }
    if (std::any_of(poly.holes_begin(), poly.holes_end(),
                    [&p](const Polygon2 &hole) {
                      return hole.bounded_side(p) == CGAL::ON_BOUNDARY;
                    })) {
      // we only want internal points
      continue;
    }
    points.push_back(v->point());
  }
  return points;
}

class VisibilityPolygonCalculator {
public:
  VisibilityPolygonCalculator(Polygon2WithHoles &poly) {
    if (!poly.outer_boundary().is_simple()) {
      throw std::runtime_error("Polygon is not simple");
    }
    polygon = poly;
    std::vector<Segment2> segments;
    if (poly.outer_boundary().area() < 0) {
      throw std::runtime_error("Polygon is not counterclockwise oriented");
    }
    for (const auto e : poly.outer_boundary().edges()) {
      auto s = e.source();
      auto t = e.target();
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
    for (const auto &hole : poly.holes()) {
      if (hole.area() > 0) {
        throw std::runtime_error("Hole is not clockwise oriented");
      }
      for (const auto e : hole.edges()) {
        auto s = e.source();
        auto t = e.target();
        auto seg = Segment2(s, t);
        segments.push_back(seg);
      }
    }
    CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());
    pl = PointLocation(env);
    auto face = env.unbounded_face();
    if (face->number_of_holes() != 1 || !face->is_unbounded()) {
      throw std::runtime_error("Bad arrangement. Could not determine polygon face.");
    }
    auto hole_it = face->holes_begin();
    assert(hole_it != face->holes_end());
    auto f = (*hole_it)->twin()->face();
    if (f->is_unbounded()) {
      throw std::runtime_error("Bad arrangement. Face should not be unbounded.");
    }
    interior_face = f;
  }

  bool is_feasible_query_point(const Point &query_point) {
    if (polygon.outer_boundary().bounded_side(query_point) ==
        CGAL::ON_UNBOUNDED_SIDE) {
      return false;
    }
    for (const auto &hole : polygon.holes()) {
      if (hole.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE) {
        return false;
      }
    }
    return true;
  }

  Face_handle get_visibility_face_handle(const Point &query_point,
                                         Arrangement_2 &output_arr) {
    /**
     * @brief The following code is slighly more complex and inefficent than
     * it may seem necessary. The reason is that it seems insufficient to just
     * provide the face handle if the query point is on the boundary of the
     * polygon. Instead, a point location has to be performed on the
     * arrangement. Still not 100% sure if there is not some bug.
     * The code looks complex because the iterators of the arrangement are
     * old school C++ iterators, which are not very nice to use. Maybe there
     * are better ways: The documentation of CGAL is old, too.
     */
    auto location = pl.locate(query_point);
    const Arrangement_2::Vertex_const_handle *v;
    const Arrangement_2::Halfedge_const_handle *e;
    const Arrangement_2::Face_const_handle *f;
    TEV tev(env);
    if (f = boost::get<Arrangement_2::Face_const_handle>(&location)) {
      if (*f != interior_face) {
        throw std::runtime_error(
            "Query point is in a face, but not in the interior face");
      }
      return tev.compute_visibility(query_point, interior_face, output_arr);
    } else if (v = boost::get<Arrangement_2::Vertex_const_handle>(&location)) {
      auto e_ = (*v)->incident_halfedges();
      if (e_->face() == interior_face) {
        return tev.compute_visibility(query_point, e_, output_arr);
      }
      ++e_;
      while (e_ != (*v)->incident_halfedges()) {
        if (e_->face() == interior_face) {
          return tev.compute_visibility(query_point, e_, output_arr);
        }
        ++e_;
      }
      throw std::runtime_error(
          "Query point is on a vertex, but not in interior face");
    } else if (e = boost::get<Arrangement_2::Halfedge_const_handle>(
                   &location)) {
      if ((*e)->face() == interior_face) {
        return tev.compute_visibility(query_point, *e, output_arr);
      } else if ((*e)->twin()->face() == interior_face) {
        return tev.compute_visibility(query_point, (*e)->twin(), output_arr);
      } else {
        throw std::runtime_error(
            "Query point is on an edge, but not the interior face");
      }
    }
    throw std::runtime_error("Unknown query point location");
  }

  Polygon2 compute_visibility_polygon(const Point &query_point) {
    if (!is_feasible_query_point(query_point)) {
      throw std::runtime_error("Query point not feasible");
    }
    Arrangement_2 output_arr;
    Face_handle fh = get_visibility_face_handle(query_point, output_arr);
    Polygon2 poly;
    std::vector<Point> points;
    auto e = fh->outer_ccb();
    points.push_back(e->source()->point());
    ++e;
    while (e != fh->outer_ccb()) {
      points.push_back(e->source()->point());
      ++e;
    }
    return Polygon2(points.begin(), points.end());
  }

  Polygon2WithHoles polygon;
  Arrangement_2 env;
  Face_handle interior_face;
  PointLocation pl;
};

template <typename EdgeHandle>
Polygon2 _boundary_to_polygon(const EdgeHandle &e) {
  Polygon2 poly;
  std::vector<Point> points;
  auto e_ = e;
  points.push_back(e_->source()->point());
  ++e_;
  while (e_ != e) {
    points.push_back(e_->source()->point());
    ++e_;
  }
  return Polygon2(points.begin(), points.end());
}

template <typename FaceHandle> Polygon2 _face_to_polygon(const FaceHandle &fh) {
  assert(!fh->is_unbounded());
  return _boundary_to_polygon(fh->outer_ccb());
}

Arrangement_2 _polygon_to_arrangement(const Polygon2WithHoles &poly) {
  // Create a new arrangement from a polygon with holes
  Arrangement_2 env;
  std::vector<Segment2> segments;
  // Add the outer boundary
  for (const auto e : poly.outer_boundary().edges()) {
    auto s = e.source();
    auto t = e.target();
    auto seg = Segment2(s, t);
    segments.push_back(seg);
  }
  // Add the holes
  for (const auto &hole : poly.holes()) {
    for (const auto e : hole.edges()) {
      auto s = e.source();
      auto t = e.target();
      auto seg = Segment2(s, t);
      segments.push_back(seg);
    }
  }
  // Insert the segments into the arrangement
  CGAL::insert(env, segments.begin(), segments.end());
  return env;
}

std::vector<Polygon2WithHoles> repair(const Polygon2WithHoles &poly) {
  // Repair a polygon with holes that is self intersecting.
  // Use arrangements to separate the polygons.

  // Create an arrangement
  Arrangement_2 env = _polygon_to_arrangement(poly);
  std::vector<Polygon2WithHoles> result;
  // Get the faces
  for (auto f = env.faces_begin(); f != env.faces_end(); ++f) {
    // Face is a polygon if it is adjacent to the unbounded face
    if (f->is_unbounded()) {
      continue;
    }
    if (!f->outer_ccb()->twin()->face()->is_unbounded()) {
      continue;
    }
    // face to polygon with holes
    // outer boundary
    auto outer_boundary = _face_to_polygon(f);
    // holes
    std::vector<Polygon2> holes;
    for (auto h = f->holes_begin(); h != f->holes_end(); ++h) {
      // h is Ccb_halfedge_circulator
      auto hole_poly = _boundary_to_polygon(*h);
      hole_poly.reverse_orientation();
      holes.push_back(hole_poly);
      assert(holes.back().area() < 0);
    }
    result.push_back(
        Polygon2WithHoles(outer_boundary, holes.begin(), holes.end()));
  }
  return result;
}

// Getting this name right is important! It has to equal the name in the
// CMakeLists.txt.
PYBIND11_MODULE(_cgal_bindings, m) {
  namespace py = pybind11;
  m.doc() = "Example of PyBind11 and CGAL."; // optional module docstring

  // Exact numbers
  py::class_<Kernel::FT>(m, "FieldNumber",
                         "A container for exact numbers in CGAL.")
      .def(py::init<long>())
      .def(py::init<double>())
      .def(py::self / Kernel::FT())
      .def(py::self + Kernel::FT())
      .def(py::self - Kernel::FT())
      .def(py::self * Kernel::FT())
      .def(py::self == Kernel::FT())
      .def(py::self < Kernel::FT())
      .def(py::self > Kernel::FT())
      .def(py::self <= Kernel::FT())
      .def(py::self >= Kernel::FT())
      .def("__float__", &CGAL::to_double<Kernel::FT>)
      .def("__str__", [](const Kernel::FT &x) {
        return std::to_string(CGAL::to_double(x));
      });

  // Points
  py::class_<Point>(m, "Point", "A 2-dimensional point")
      .def(py::init<long, long>())
      .def(py::init<double, double>())
      .def(py::init<Kernel::FT, Kernel::FT>())
      .def("x", [](const Point &p) { return p.x(); })
      .def("y", [](const Point &p) { return p.y(); })
      .def("__len__", [](const Point &self) { return 2; })
      .def("__item__",
           [](const Point &self, int i) {
             if (i == 0) {
               return self.x();
             } else if (i == 1) {
               return self.y();
             }
             throw std::out_of_range("Only 0=x and 1=y.");
           })
      .def(py::self == Point())
      .def("__str__", [](const Point &p) {
        return fmt::format("({}, {})", CGAL::to_double(p.x()),
                           CGAL::to_double(p.y()));
      });

  // Polygons
  py::class_<Polygon2>(m, "Polygon", "A simple polygon in CGAL.")
      .def(py::init<>())
      .def(py::init([](const std::vector<Point> &vertices) {
        return std::make_unique<Polygon2>(vertices.begin(), vertices.end());
      }))
      .def("boundary",
           [](const Polygon2 &poly) {
             std::vector<Point> points;
             std::copy(poly.begin(), poly.end(), std::back_inserter(points));
             return points;
           })
      .def("is_simple", &Polygon2::is_simple)
      .def("contains",
           [](const Polygon2 &self, const Point &p) {
             return self.bounded_side(p) != CGAL::ON_UNBOUNDED_SIDE;
           })
      .def("on_boundary",
           [](const Polygon2 &self, const Point &p) {
             return self.bounded_side(p) == CGAL::ON_BOUNDARY;
           })
      .def("interior_sample_points",
           [](Polygon2 &self) { return get_sample_points(self); })
      .def("area", [](const Polygon2 &poly) { return poly.area(); });

  py::class_<Polygon2WithHoles>(m, "PolygonWithHoles",
                                "A polygon with holes in CGAL.")
      .def(py::init([](const Polygon2 &outer,
                       const std::vector<Polygon2> &holes) {
        for (const auto &hole_poly : holes) {
          if (hole_poly.area() >= 0) {
            throw std::runtime_error("Hole is not clockwise oriented");
          }
        }
        if (outer.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        return Polygon2WithHoles(outer, holes.begin(), holes.end());
      }))
      .def(py::init([](const Polygon2 &outer) {
        if (outer.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        return Polygon2WithHoles(outer);
      }))
      .def(py::init([](const std::vector<Point> &outer_vertices) {
        auto poly = Polygon2(outer_vertices.begin(), outer_vertices.end());
        return Polygon2WithHoles(poly);
      }))
      .def(py::init([](const std::vector<Point> &outer_vertices,
                       const std::vector<std::vector<Point>> &hole_vertices) {
        auto poly = Polygon2(outer_vertices.begin(), outer_vertices.end());
        if (poly.area() <= 0) {
          throw std::runtime_error("Polygon is not counterclockwise oriented");
        }
        std::vector<Polygon2> holes;
        for (auto hole_boundary : hole_vertices) {
          auto hole_poly = Polygon2(hole_boundary.begin(), hole_boundary.end());
          if (hole_poly.area() >= 0) {
            throw std::runtime_error("Hole is not clockwise oriented");
          }
          holes.push_back(hole_poly);
        }
        return Polygon2WithHoles(poly, holes.begin(), holes.end());
      }))
      .def(py::self == Polygon2WithHoles())
      .def(
          "outer_boundary",
          [](const Polygon2WithHoles &self) { return self.outer_boundary(); },
          "Returns a list with all holes (simple polygons).")
      .def("holes",
           [](const Polygon2WithHoles &poly) {
             // Copy the holes into a vector so that PyBind11 can return it
             // as a Python list.
             std::vector<Polygon2> holes;
             std::copy(poly.holes_begin(), poly.holes_end(),
                       std::back_inserter(holes));
             return holes;
           })
      .def("contains",
           [](Polygon2WithHoles &self, const Point &p) {
             if (self.outer_boundary().bounded_side(p) ==
                 CGAL::ON_UNBOUNDED_SIDE) {
               return false;
             }
             for (auto hole : self.holes()) {
               if (hole.bounded_side(p) == CGAL::ON_BOUNDED_SIDE) {
                 return false;
               }
             }
             return true;
           })
      .def("on_boundary",
           [](Polygon2WithHoles &self, const Point &p) {
             if (self.outer_boundary().bounded_side(p) == CGAL::ON_BOUNDARY) {
               return true;
             }
             for (const auto &hole : self.holes()) {
               if (hole.bounded_side(p) == CGAL::ON_BOUNDARY) {
                 return true;
               }
             }
             return false;
           })
      .def("interior_sample_points",
           [](Polygon2WithHoles &self) {
             return get_sample_points_with_holes(self);
           })
      .def("area",
           [](const Polygon2WithHoles &poly) {
             auto area = poly.outer_boundary().area();
             for (const auto &hole : poly.holes()) {
               area -= hole.area();
             }
             return CGAL::to_double(area);
           })
      .def(
          "join",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            Polygon2WithHoles joined;
            if (CGAL::join(self, other, joined,
                           /*UsePolylines=*/CGAL::Tag_false{})) {
              result.push_back(joined);
            } else {
              result.push_back(self);
              result.push_back(other);
            }
            return result;
          },
          "Joins two polygons (with holes) into a list of polygons (with "
          "holes).")
      .def(
          "intersection",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            CGAL::intersection(self, other, std::back_inserter(result),
                               /*UsePolylines=*/CGAL::Tag_false{});
            return result;
          },
          "Computes the intersection of two polygons (with holes). Returns a "
          "list of polygons (with holes).")
      .def(
        "do_intersect",
        [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
          return CGAL::do_intersect(self, other);
        },
        "Check if two polygons (with holes) intersect. Returns a boolean.")
      .def(
          "difference",
          [](const Polygon2WithHoles &self, const Polygon2WithHoles &other) {
            std::vector<Polygon2WithHoles> result;
            CGAL::difference(self, other, std::back_inserter(result),
                             /*UsePolylines=*/CGAL::Tag_false{});
            return result;
          },
          "Removes the area of the other polygon. Returns a list of polygons "
          "(with holes).");

  py::class_<VisibilityPolygonCalculator>(
      m, "VisibilityPolygonCalculator",
      "A class to compute visibility polygons.")
      .def(py::init<Polygon2WithHoles &>())
      .def("compute_visibility_polygon",
           &VisibilityPolygonCalculator::compute_visibility_polygon,
           py::arg("query_point"),
           "Compute the visibility polygon for a query point.")
      .def("is_feasible_query_point",
           &VisibilityPolygonCalculator::is_feasible_query_point,
           py::arg("query_point"),
           "Check if the query point is within the polygon.");

  m.def("repair", &repair,
        "Repair a polygon with holes that is self "
        "intersecting. Returns a list of polygons with "
        "holes.");
}
