#pragma once
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// This is necessary to build on appveyor for some reason
#define CGAL_EIGEN3_ENABLED 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>


namespace py = pybind11;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_2                 Point2;
typedef Kernel::Point_3                 Point3;
typedef Kernel::Vector_2                Vector2;
typedef Kernel::Vector_3                Vector3;

typedef CGAL::Surface_mesh<Point3>      Mesh3;
typedef Mesh3::Vertex_index             V;
typedef Mesh3::Face_index               F;
typedef Mesh3::Halfedge_index           H;
typedef Mesh3::Edge_index               E;
