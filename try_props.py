from seagullmesh import Mesh3
import numpy as np
from numpy import array, cos, sin
from numpy.testing import assert_array_equal


def tetrahedron(scale=1.0, rot_z=0.0):
    verts = scale * array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype='float')

    if rot_z:
        rot = array([[cos(rot_z), -sin(-rot_z), 0], [sin(rot_z), cos(-rot_z), 0], [0, 0, 1]])
        verts = verts @ rot.T

    faces = array([[2, 1, 0], [2, 3, 1], [3, 2, 0], [1, 3, 0]], dtype='int')
    return verts, faces


mesh = Mesh3.from_polygon_soup(*tetrahedron())
# mesh.face_data.add_property('what', default=np.ulonglong(2))
mesh.face_data.add_property('what', default=np.int_(2))
# assert_array_equal(mesh.face_data['what'][:], [-1, -1, -1, -1])
print(mesh.face_data['what'].pmap)

# mesh.face_data.add_property('huh', default=np.uint(2))
# assert_array_equal(mesh.face_data['huh'][:], [2, 2, 2, 2])
# print(mesh.face_data['huh'].pmap)
