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


def test_from_polygon_soup():
    verts, faces = tetrahedron()
    mesh = Mesh3.from_polygon_soup(verts, faces)
    assert mesh.n_vertices == 4 and mesh.n_faces == 4


@pytest.mark.parametrize('file', ['armadillo.off', 'sphere.ply'])
def test_from_file(file):
    file = Path(__file__).parent / 'assets' / file
    assert file.exists()
    _mesh = Mesh3.from_file(str(file))


@pytest.mark.parametrize('ext', ['ply', 'off'])
def test_to_file(ext):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    with TemporaryDirectory() as d:
        file = str(Path(d) / f'mesh.{ext}')
        mesh.to_file(file)


def test_pyvista_roundtrip():
    from pyvista import Sphere
    pvmesh0 = Sphere().clean().triangulate()
    mesh = Mesh3.from_pyvista(pvmesh0)
    _pvmesh1 = mesh.to_pyvista()


@pytest.mark.parametrize(
    ['cls', 'default'],
    [
        (props.VertBoolPropertyMap, False),
        (props.VertIntPropertyMap, 0),
        (props.VertDoublePropertyMap, 0.0),
        (props.VertPoint2PropertyMap, Point2(0, 0)),
    ]
)
def test_property_map_type_casting(cls, default):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    d = mesh.vertex_data
    pmap = d.add_property('foo', default=default)
    assert isinstance(pmap.pmap, cls)


@pytest.mark.parametrize('key_type', ['vertex', 'face', 'edge', 'halfedge'])
@pytest.mark.parametrize('val_type', [int, bool, float])
def test_scalar_properties(key_type, val_type):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    d: MeshData = getattr(mesh, key_type + '_data')

    d['foo'] = full(d.n_mesh_keys, val_type(0))

    keys = d.mesh_keys
    key = keys[0]
    d['foo'][key] = val_type(1)

    assert 'foo' in d
    assert 'bar' not in d

    val = d['foo'][key]
    assert val == val_type(1)

    d['foo'][keys[:2]] = [val_type(1), val_type(1)]
    assert d['foo'][keys[0]] == val_type(1) and d['foo'][keys[1]] == val_type(1)

    d.remove_property('foo')
    assert 'foo' not in d.mesh_keys
    assert 'foo' not in d


@pytest.mark.parametrize('key_type', ['vertex', 'face', 'edge', 'halfedge'])
@pytest.mark.parametrize('val_type', [Point2, Point3, Vector2, Vector3])
def test_array_properties(key_type, val_type):
    mesh = Mesh3.from_polygon_soup(*tetrahedron())
    d: MeshData = getattr(mesh, key_type + '_data')

    ndims = int(val_type.__name__[-1])
    default = val_type(*[0.0 for _ in range(ndims)])
    d.add_property('foo', default=default)

    nkeys = d.n_mesh_keys
    data = np.random.uniform(-1, 1, (nkeys, ndims))
    d['foo'] = data

    assert_array_equal(data, d['foo'][:])


def corefine_meshes():
    m1 = Mesh3.from_polygon_soup(*tetrahedron())
    m2 = Mesh3.from_polygon_soup(*tetrahedron(scale=0.9, rot_z=pi/3))
    return m1, m2


def test_corefine():
    m1, m2 = corefine_meshes()
    nv1_orig, nv2_orig = m1.n_vertices, m2.n_vertices
    m1.corefine(m2)

    nv1, nv2 = m1.n_vertices, m2.n_vertices
    assert nv1 > nv1_orig and nv2 > nv2_orig


@pytest.mark.parametrize('op', ['union', 'intersection', 'difference'])
@pytest.mark.parametrize('inplace', [False, True])
def test_boolean_ops(op, inplace):
    m1, m2 = corefine_meshes()
    _m3 = getattr(m1, op)(m2, inplace=inplace)


@pytest.fixture()
def armadillo():
    file = Path(__file__).parent / 'assets' / 'armadillo.off'
    assert file.exists()
    return Mesh3.from_file(str(file))


def test_estimate_geodesic_distance_source_vert(armadillo):
    armadillo.estimate_geodesic_distances(armadillo.vertices[0], 'distances')
    assert (armadillo.vertex_data['distances'] > 0).any()


def test_estimate_geodesic_distance_source_verts(armadillo):
    armadillo.estimate_geodesic_distances(armadillo.vertices[:3], 'distances')
    assert (armadillo.vertex_data['distances'] > 0).any()
