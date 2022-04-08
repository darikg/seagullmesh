#include "util.hpp"

std::vector<Point3> array_to_points_3(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 3) {
        throw std::runtime_error("vertices need to be 3 dimensional");
    }
    const size_t nv = v.shape(0);
    std::vector<Point3> points;
    points.reserve(nv);
    for (size_t i = 0; i < nv; i++) {
        points.emplace_back(Point3(v(i, 0), v(i, 1), v(i, 2)));
    }
    return points;
}

std::vector<Point2> array_to_points_2(const py::array_t<double> &verts) {
    auto v = verts.unchecked<2>();
    if (v.shape(1) != 2) {
        throw std::runtime_error("vertices need to be 2 dimensional");
    }
    const size_t nv = v.shape(0);
    std::vector<Point2> points;
    points.reserve(nv);
    for (size_t i = 0; i < nv; i++) {
        points.emplace_back(Point2(v(i, 0), v(i, 1)));
    }
    return points;
}

template<typename Point>
py::array_t<double, py::array::c_style> points_d_to_array(const std::vector<Point>& points, const size_t ndims) {
    const size_t np = points.size();
    py::array_t<double, py::array::c_style> points_out({np, ndims});
    auto r = points_out.mutable_unchecked<2>();
    for (auto i = 0; i < np; i++) {
        for (auto j = 0; j < ndims; j++) {
            r(i, j) = CGAL::to_double(points[i][j]);
        }
    }
    return points_out;
}

py::array_t<double, py::array::c_style> points_to_array(const std::vector<Point3>& points) {
    return points_d_to_array<Point3>(points, 3);
}

py::array_t<double, py::array::c_style> points_to_array(const std::vector<Point2>& points) {
    return points_d_to_array<Point2>(points, 2);
}
