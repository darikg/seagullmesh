import pytest
import os

from tempfile import TemporaryDirectory, TemporaryFile
from pathlib import Path

import numpy as np
from numpy import array, cos, sin,  pi, ones, arange, uint64, full, zeros
from numpy.testing import assert_array_equal

from seagullmesh import sgm, Mesh3, MeshData, Point2, Point3, Vector2, Vector3
props = sgm.properties


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot.T

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')
    return verts, faces


def tetrahedron_mesh() -> Mesh3:
    return Mesh3.from_polygon_soup(*tetrahedron())


# @pytest.mark.parametrize('key_type', ['vertex', 'face', 'edge', 'halfedge'])
# @pytest.mark.parametrize('val_type', [int, bool, float])
# def test_scalar_properties(key_type, val_type):
#     mesh = Mesh3.from_polygon_soup(*tetrahedron())
#     data: MeshData = getattr(mesh, key_type + '_data')
#
#     data['foo'] = full(data.n_mesh_keys, val_type(0))
#
#     keys = data.mesh_keys
#     key = keys[0]
#     data['foo'][key] = val_type(1)
#
#     assert 'foo' in data
#     assert 'bar' not in data
#
#     val = data['foo'][key]
#     assert val == val_type(1)
#
#     data['foo'][keys[:2]] = [val_type(1), val_type(1)]
#     assert data['foo'][keys[0]] == val_type(1) and data['foo'][keys[1]] == val_type(1)
#
#     data.remove_property('foo')
#     assert 'foo' not in data.mesh_keys
#     assert 'foo' not in data


def test_half_edge():
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    # halfedges = mesh._mesh.halfedges
    # print(list(mesh._mesh.first_halfedge))
    print(sgm.mesh.Halfedge(0))
    # pmap = sgm.properties.HalfedgeBoolPropertyMap(mesh.mesh, 'foo', False)
    # pmap[halfedges[2]] = True
    # assert pmap[halfedges[2]]
