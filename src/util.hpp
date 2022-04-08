#pragma once
#include "seagullmesh.hpp"

std::vector<Point3> array_to_points_3(const py::array_t<double> &verts);
std::vector<Point2> array_to_points_2(const py::array_t<double> &verts);

py::array_t<double, py::array::c_style> points_to_array(const std::vector<Point3>& points);
py::array_t<double, py::array::c_style> points_to_array(const std::vector<Point2>& points);
