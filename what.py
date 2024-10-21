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


mesh = Mesh3.from_polygon_soup(*tetrahedron())
pmap = mesh.halfedge_data.add_property('foo', default=False)
print(sgm.mesh.Halfedge(0))
# print(mesh.mesh.first_halfedge)
print(mesh.mesh.halfedges)
print(pmap.pmap[mesh.mesh.halfedges])
