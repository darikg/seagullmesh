from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import contextmanager
from functools import cached_property
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING, Union, Sequence, TypeVar, overload, Tuple, \
    Generic, Iterator, Type, Dict, Literal, Callable, Sized

import numpy as np
from seagullmesh._seagullmesh.mesh import (
    Mesh3 as _Mesh3,
    Point2, Point3, Vector2, Vector3,
    Vertex, Face, Edge, Halfedge,
)
from typing_extensions import Self

from seagullmesh import _seagullmesh as sgm
from ._version import version_info, __version__  # noqa

if TYPE_CHECKING:
    try:
        import pyvista as pv  # noqa
    except ImportError:
        pv = None

_IndexTypes = Vertex | Face | Edge | Halfedge
_CppIndicesTypes = sgm.mesh.Vertices | sgm.mesh.Faces | sgm.mesh.Edges | sgm.mesh.Halfedges
TIndex = TypeVar('TIndex', Vertex, Face, Edge, Halfedge)
TIndices = TypeVar('TIndices', sgm.mesh.Vertices, sgm.mesh.Faces, sgm.mesh.Edges, sgm.mesh.Halfedges)


class Indices(Generic[TIndex, TIndices], Sequence[TIndex], Sized):
    index_type: Type[TIndex]  # Set by subclass
    indices_type: Type[_CppIndicesTypes]  # The C++ indices class, set by subclass

    # Updated in __init__subclass__
    _cpp_indices_to_py_indices: Dict[Type[_CppIndicesTypes], Type[_PyIndicesTypes]] = {}

    def __init__(self, mesh: Mesh3, indices: TIndices | Sequence[TIndex]):
        self.mesh = mesh

        if isinstance(indices, _CppIndicesTypes):
            if isinstance(indices, self.indices_type):
                self.indices = indices
            else:
                msg = f"Can't construct {type(self).__name__} from {type(indices).__name__}"
                raise TypeError(msg)
        else:
            # A sequence e.g. list[Face] we can coerce into Indices<Face>?
            # Should raise a pybind11 type error otherwise
            self.indices = self.indices_type(indices)

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__()
        Indices._cpp_indices_to_py_indices[cls.indices_type] = cls  # type: ignore

    @classmethod
    @abstractmethod
    def all_indices(cls, mesh: Mesh3) -> Self: ...

    @classmethod
    @abstractmethod
    def n_mesh_keys(cls, mesh: Mesh3) -> int: ...

    @staticmethod
    def py_indices_type(cpp_indices_type: _CppIndicesTypes) -> Type[_PyIndicesTypes]:
        try:
            return Indices._cpp_indices_to_py_indices[cpp_indices_type]
        except KeyError:
            raise TypeError(f"{cpp_indices_type} is not a seagullmesh C++ indices type")

    @staticmethod
    def from_indices(mesh: Mesh3, indices: TIndices) -> Indices:
        return Indices._cpp_indices_to_py_indices[type(indices)](mesh, indices)

    @property
    def _array(self) -> np.ndarray:  # array of ints
        return self.indices.indices  # mapped to Indices.get_indices() on the C++ side

    def _with_array(self, arr: np.NDArray[np.uint32]) -> Self:
        indices = self.indices_type(arr)
        return type(self)(self.mesh, indices)

    def __repr__(self) -> str:
        return f'{self.indices_type.__name__}(n={len(self)})'

    def __len__(self) -> int:
        return len(self._array)

    def __eq__(self, other: TIndex | Self) -> np.ndarray:
        # e.g. (these_faces == that_face) -> array[bool]
        if isinstance(other, self.index_type):
            return self._array == other.to_int()

        # e.g. (these_faces == those_faces) -> array[bool]
        elif isinstance(other, Indices) and (other.index_type is self.index_type):
            return self._array == other._array
        else:
            msg = f"Can only compare indices of the same type, got self={self} and other={other}"
            raise TypeError(msg)

    def __ne__(self, other: TIndex | Self) -> np.ndarray:
        return np.logical_not((self == other))

    def __iter__(self) -> Iterator[TIndex]:
        for i in self._array:
            yield self.index_type(i)

    def copy(self) -> Self:
        return self._with_array(self._array.copy())

    @overload
    def __getitem__(self, item: int | np.integer) -> TIndex: ...  # Vertex, Face, Edge, Halfedge

    @overload
    def __getitem__(self, item: slice | Sequence[int] | np.ndarray) -> Self: ...  # slice self to get another Self

    def __getitem__(self, item):
        if isinstance(item, (int, np.integer)):
            return self.index_type(self._array[item])  # Convert int to descriptor
        else:
            # As long as it can slice an array we're happy
            return self._with_array(self._array[item])

    @overload
    def __setitem__(self, item: int, value: TIndex): ...  # self.idxs[8] = some_vertex

    @overload
    def __setitem__(self, item: Sequence, value: Self): ...  # self.idxs[1:]] = self.idxs[1:][::-1]

    def __setitem__(self, item, value):
        if isinstance(item, int) and isinstance(value, self.index_type):
            # faces[0] = some_face
            self._array[item] = value.to_int()
            return

        if isinstance(value, type(self)):
            if value.mesh is self.mesh:
                self._array[item] = value._array
                return
            else:
                raise ValueError("Assigning indices between different meshes, this is probably a mistake")

        msg = f"Can only assign indices of the same type, got self={self} and value={value}"
        raise TypeError(msg)

    def unique(self) -> Self:
        return self._with_array(np.unique(self._array))


class Vertices(Indices[Vertex, sgm.mesh.Vertices]):
    index_type = Vertex
    indices_type = sgm.mesh.Vertices

    @classmethod
    def all_indices(cls, mesh: Mesh3) -> Self:
        return cls(mesh, mesh.mesh.vertices)

    @classmethod
    def n_mesh_keys(cls, mesh: Mesh3) -> int:
        return mesh.mesh.n_vertices

    def adjacent_faces(self) -> Faces:
        return Faces(self.mesh, sgm.connected.vertices_to_faces(self.mesh.mesh, self.indices))

    def adjacent_edges(self) -> Edges:
        return Edges(self.mesh, sgm.connected.edges_to_faces(self.mesh.mesh, self.indices))

    def degrees(self) -> np.ndarray:
        return sgm.connected.vertex_degrees(self.mesh, self.indices)

    def points(self) -> np.ndarray:
        return self.mesh.vertex_point_map[self]

    def normals(self) -> np.ndarray:
        return sgm.geometry.vertex_normals(self.mesh.mesh, self.indices)


class Faces(Indices[Face, sgm.mesh.Faces]):
    index_type = Face
    indices_type = sgm.mesh.Faces

    @classmethod
    def all_indices(cls, mesh: Mesh3) -> Self:
        return cls(mesh, mesh.mesh.faces)

    @classmethod
    def n_mesh_keys(cls, mesh: Mesh3) -> int:
        return mesh.mesh.n_faces

    def degrees(self) -> np.ndarray:
        return sgm.connected.face_degrees(self.mesh, self.indices)

    def construct_points(
            self,
            bary_coords: np.ndarray,
            vert_points: str | PropertyMap[Vertex, Point2 | Point3] | None = None,
    ) -> np.ndarray:
        """Construct a set of points from face barycentric coordinates

        `bary_coords` must be of shape (len(faces), 3)

        By default, points are constructed using the default mesh vertex points. An optional vertex point map of value
        Point2 or Point3 can also be supplied. The returned array if of shape (len(faces), 2 or 3) as appropriate.
        """
        vpm = self.mesh.get_vertex_point_map(vert_points)
        return sgm.locate.construct_points(self.mesh.mesh, self.indices, bary_coords, vpm.pmap)

    def triangle_soup(self) -> np.ndarray:
        return sgm.io.triangle_soup(self.mesh.mesh, self.indices)  # TODO index, index_map

    def normals(self) -> np.ndarray:
        """(nf, 3) array of face normal vectors"""
        return sgm.geometry.face_normals(self.mesh.mesh, self.indices)

    def areas(self) -> np.ndarray:
        """(nf,) array of face areas"""
        return sgm.geometry.face_areas(self.mesh.mesh, self.indices)

    def adjacent_edges(self) -> Edges:
        return Edges(self.mesh, sgm.connected.faces_to_edges(self.mesh.mesh, self.indices))

    def adjacent_vertices(self) -> Vertices:
        return Vertices(self.mesh, sgm.connected.faces_to_vertices(self.mesh.mesh, self.indices))

    def is_null(self) -> np.ndarray:
        return self == Mesh3.null_face

    def contains_point(self, point: Point3) -> np.ndarray:
        return sgm.locate.is_point_in_faces(self.mesh.mesh, point, self.indices)


class Edges(Indices[Edge, sgm.mesh.Edges]):
    index_type = Edge
    indices_type = sgm.mesh.Edges

    @classmethod
    def all_indices(cls, mesh: Mesh3) -> Self:
        return cls(mesh, mesh.mesh.edges)

    @classmethod
    def n_mesh_keys(cls, mesh: Mesh3) -> int:
        return mesh.mesh.n_edges

    def lengths(self) -> np.ndarray:
        return sgm.geometry.edge_lengths(self.mesh.mesh, self.indices)


class Halfedges(Indices[Halfedge, sgm.mesh.Halfedges]):
    index_type = Halfedge
    indices_type = sgm.mesh.Halfedges

    @classmethod
    def all_indices(cls, mesh: Mesh3) -> Self:
        return cls(mesh, mesh.mesh.halfedges)

    @classmethod
    def n_mesh_keys(cls, mesh: Mesh3) -> int:
        return mesh.mesh.n_halfedges


_PyIndicesTypes = Vertices | Faces | Edges | Halfedges


class Mesh3:
    null_vertex: Vertex = _Mesh3.null_vertex
    null_face: Face = _Mesh3.null_face
    null_edge: Edge = _Mesh3.null_edge
    # null_halfedge: Halfedge = _Mesh3.null_halfedge

    def __init__(self, mesh: _Mesh3 | None = None):
        """Construct a python-wrapped mesh"""

        """The C++ mesh object being wrapped"""
        self.mesh = mesh if mesh else _Mesh3()

        if hasattr(sgm, 'properties'):
            # Allow installing seagullmesh without the properties modules
            self.vertex_data = self.vd = VertexData(self)
            self.face_data = self.fd = FaceData(self)
            self.edge_data = self.ed = EdgeData(self)
            self.halfedge_data = self.hd = HalfedgeData(self)

    n_vertices = nv = property(lambda self: self.mesh.n_vertices)
    n_faces = nf = property(lambda self: self.mesh.n_faces)
    n_edges = ne = property(lambda self: self.mesh.n_edges)
    n_halfedges = nh = property(lambda self: self.mesh.n_halfedges)

    vertices = vs = property(lambda self: Vertices.all_indices(self))
    faces = fs = property(lambda self: Faces.all_indices(self))
    edges = es = property(lambda self: Edges.all_indices(self))
    halfedges = hs = property(lambda self: Halfedges.all_indices(self))

    is_valid = property(lambda self: self.mesh.is_valid(False))

    def is_mesh_valid(self, verbose=False):
        return self.mesh.mesh(verbose)

    @property
    def is_closed(self):
        return self.mesh.is_closed()

    @property
    def has_garbage(self) -> bool:
        return self.mesh.has_garbage

    def collect_garbage(self) -> None:
        self.mesh.collect_garbage()

    @cached_property
    def vertex_point_map(self) -> PropertyMap[Vertex, Point3]:
        # This is maintained internally by cgal
        return self.vertex_data.find(sgm.properties.V_Point3_PropertyMap, name='v:point')

    @cached_property
    def _vertex_index_map(self) -> VertexIndexMap:
        cpp_pmap = sgm.properties.V_uint32_PropertyMap.get_or_create(self.mesh, 'v:index', 0)
        return self.vertex_data.wrap(cpp_pmap, dtype_name='uint32')

    @property
    def vertex_index_map(self) -> VertexIndexMap:
        pmap = self._vertex_index_map
        pmap.pmap.index(self.mesh)  # Make sure indices are up-to-date
        return pmap

    def iter_meshdata(self) -> Iterator[MeshData]:
        yield self.vertex_data
        yield self.face_data
        yield self.edge_data
        yield self.halfedge_data

    def copy(self) -> Mesh3:
        """Deep-copy the mesh and all its properties"""
        out = Mesh3(sgm.mesh.Mesh3(self.mesh))

        # The properties have been copied, just need to create new pmap references
        for src_data, dest_data in zip(self.iter_meshdata(), out.iter_meshdata()):
            for name, src_wrapper in src_data.items():
                wrapped = dest_data.find(
                    pmap_cls=type(src_wrapper.pmap),
                    name=name,
                    wrapper_cls=type(src_wrapper),
                    dtype_name=src_wrapper.dtype_name,
                )
                assert isinstance(wrapped, PropertyMap)
                dest_data[name] = wrapped

        return out

    def add(self, other: Mesh3, check_properties=False, inplace=False) -> Mesh3:
        if check_properties:
            for d_self, d_other in zip(self.iter_meshdata(), other.iter_meshdata()):
                d_self.check_has_same_properties(d_other)

        out = self if inplace else self.copy()
        out.mesh += other.mesh
        return out

    @staticmethod
    def from_polygon_soup(verts: np.ndarray, faces: np.ndarray, orient=True, validate=False) -> Mesh3:
        """Constructs a surface mesh from vertices (nv * 3) and faces (nf * 3) arrays

        If `orient` is True (default), the faces are reindexed to represent a consistent manifold surface.
        """
        mesh = sgm.io.polygon_soup_to_mesh3(verts, faces, orient, validate)
        return Mesh3(mesh)

    def to_polygon_soup(self) -> Tuple[np.ndarray, np.ndarray]:
        """Returns vertices (nv * 3) and faces (nf * 3) array"""
        return sgm.io.mesh3_to_polygon_soup(self.mesh)

    def triangle_soup(self) -> np.ndarray:
        return sgm.io.triangle_soup(self.mesh, self.vertex_index_map.pmap)

    def to_edge_soup(
            self,
            vim: VertexIndexMap | None = None,
            edges_only: bool = False,
    ) -> np.ndarray | Tuple[np.ndarray, np.ndarray]:
        vim = vim or self.vertex_index_map
        edges = sgm.io.edge_soup(self.mesh, vim.pmap)
        if edges_only:
            return edges

        pts = sgm.io.point_soup(self.mesh)
        return pts, edges

    def point_soup(self) -> np.ndarray:
        return sgm.io.point_soup(self.mesh)

    @staticmethod
    def from_file(filename: str) -> Mesh3:
        mesh = sgm.io.load_mesh_from_file(filename)
        return Mesh3(mesh)

    def to_file(self, filename: str):
        ext = Path(filename).suffix
        if ext == '.ply':
            sgm.io.write_ply(self.mesh, filename)
        elif ext == '.off':
            sgm.io.write_off(self.mesh, filename)
        else:
            raise ValueError(f"Unsupported format '{ext}'")

    @staticmethod
    def icosahedron(center: np.ndarray | Sequence[float] = (0, 0, 0), radius: float = 1.0) -> Mesh3:
        out = Mesh3()
        sgm.mesh.add_icosahedron(out.mesh, *center, radius)
        return out

    @staticmethod
    def pyramid(
            base_center: np.ndarray | Sequence[float] | Point3 = (0, 0, 0),
            n_base_pts: int = 4,
            height: float = 1.0,
            radius: float = 1.0,
            closed: bool = True,
    ):
        out = Mesh3()
        base_center = Point3(*base_center) if not isinstance(base_center, Point3) else base_center
        sgm.mesh.add_pyramid(out.mesh, n_base_pts, base_center, height, radius, closed)
        return out

    @staticmethod
    def grid(
            ni: int,
            nj: int,
            calculator: Callable[[int, int], Point3] = lambda i, j: Point3(i, j, 0),
            triangulated: bool = False,
    ):
        out = Mesh3()
        sgm.mesh.add_grid(out.mesh, ni, nj, calculator, triangulated)
        return out

    def transform(
            self,
            transform: np.ndarray,
            vertex_point_map: str | PropertyMap[Vertex, Point3] | None = None,
            inplace=True,
    ) -> Mesh3:
        out = self if inplace else self.copy()
        vpm = out.get_vertex_point_map(vertex_point_map)
        sgm.geometry.transform(out.mesh, transform, vpm.pmap)
        return out

    def scale(self, scale: float | Sequence[float], inplace=True) -> Mesh3:
        transform = np.diag(np.broadcast_to(scale, (3, )))
        return self.transform(transform, inplace=inplace)

    def volume(self) -> float:
        return sgm.geometry.volume(self.mesh)

    def area(self) -> float:
        return sgm.geometry.area(self.mesh)

    def centroid(self, vertex_point_map: str | PropertyMap[Vertex, Point3] | None = None) -> float:
        vpm = self.get_vertex_point_map(vertex_point_map)
        return sgm.geometry.centroid(self.mesh, vpm.pmap)

    def bounding_box(self) -> sgm.mesh.BoundingBox3:
        return sgm.geometry.bounding_box(self.mesh)

    @property
    def does_bound_a_volume(self) -> bool:
        return sgm.geometry.does_bound_a_volume(self.mesh)

    @property
    def is_outward_oriented(self) -> bool:
        return sgm.geometry.is_outward_oriented(self.mesh)

    def reverse_face_orientations(self, faces: Faces | None = None) -> None:
        if faces is None:
            sgm.geometry.reverse_face_orientations(self.mesh)
        else:
            sgm.geometry.reverse_face_orientations(self.mesh, faces)

    @staticmethod
    def from_pyvista(
            polydata: pv.PolyData,
            orient=True,
            vertex_data: Sequence[str] | Literal['all'] = (),
            face_data: Sequence[str] | Literal['all'] = (),
    ) -> Mesh3:
        """Converts a `pyvista.PolyData` object to a surface mesh.

        All point/cell data is ignored unless property names to copy are specified.
        """
        from vtkmodules.util.numpy_support import vtk_to_numpy
        cells = vtk_to_numpy(polydata.GetPolys().GetConnectivityArray())
        faces = cells.reshape(-1, 3)
        out = Mesh3.from_polygon_soup(polydata.points, faces, orient=orient)

        if vertex_data:
            keys = polydata.point_data.keys() if vertex_data == 'all' else vertex_data
            for k in keys:
                out.vertex_data[k] = polydata.point_data[k]

        if face_data:
            keys = polydata.cell_data.keys() if face_data == 'all' else face_data
            for k in keys:
                out.face_data[k] = polydata.cell_data[k]

        return out

    def to_pyvista(
            self,
            data: Literal[True] | None = None,
            vertex_data: Literal[True] | Sequence[str] = (),
            face_data: Literal[True] | Sequence[str] = (),
    ) -> pv.PolyData:
        """Returns the mesh as a `pyvista.PolyData` object.

        By default, vertex and cell data is ignored -- specify vertex_data and cell_data as a list of keys
        naming property maps to copy, or the string 'all' for all of them.
        """
        import pyvista as pv
        verts, faces = sgm.io.mesh3_to_polygon_soup(self.mesh)

        if len(set(len(f) for f in faces)) == 1:
            mesh = pv.PolyData.from_regular_faces(verts, faces)
        else:
            mesh = pv.PolyData.from_irregular_faces(verts, faces)

        if data:
            vertex_data = face_data = True

        if vertex_data:
            keys = self.vertex_data.keys() if vertex_data is True else vertex_data
            vertices = self.vertices
            for k in keys:
                mesh.point_data[k] = self.vertex_data[k][vertices]

        if face_data:
            keys = self.face_data.keys() if face_data is True else face_data
            faces = self.faces
            for k in keys:
                mesh.cell_data[k] = self.face_data[k][faces]

        return mesh

    def to_pyvista_edges(self, edge_data: Literal[True] | Sequence[str] = ()) -> pv.PolyData:
        import pyvista as pv
        mesh = pv.PolyData()
        mesh.points, edges = self.to_edge_soup()
        mesh.lines = pv.CellArray.from_regular_cells(edges)
        eidx = self.edges
        keys = self.edge_data.keys() if edge_data is True else edge_data
        for k in keys:
            mesh.cell_data[k] = self.edge_data[k][eidx]

        return mesh

    def plot(self, **kwargs):
        return self.to_pyvista(vertex_data=True, face_data=True).plot(**kwargs)

    def corefiner(self, other: Mesh3, **kwargs):
        from .corefine import Corefiner
        return Corefiner(mesh0=self, mesh1=other, **kwargs)

    def get_vertex_point_map(self, vpm: str | VertexPointMap | None) -> VertexPointMap:
        if vpm is None:
            return self.vertex_point_map
        else:
            return self.vertex_data.check(vpm)

    def aabb_tree(self, vpm: str | VertexPointMap | None = None):
        """Construct an axis-aligned bounding box tree for accelerated spatial queries

        By default, the AABB tree is constructed for the default mesh vertex locations,
        but also accepts a vertex property map storing Point2 or Point3 locations.
        """
        from .aabb import AabbTree
        tree = sgm.locate.aabb_tree(self.mesh, self.get_vertex_point_map(vpm).pmap)
        return AabbTree(self, tree)

    def estimate_geodesic_distances(
            self,
            src: Union[Vertex, Vertices],
            distance_prop: str | PropertyMap[Vertex, float],
    ) -> PropertyMap[Vertex, float]:
        """Estimates the geodesic distance from the source vertex/vertices to all vertices in the mesh

        Estimated distances are stored in the supplied vertex property map.
        """
        distances = self.vertex_data.get(distance_prop, default=0.0)
        self.mesh.estimate_geodesic_distances(distances.pmap, src)
        return distances

    def remesh(
            self,
            target_edge_length: float,
            n_iter: int = 1,
            protect_constraints=False,
            vertex_constrained: str | PropertyMap[Vertex, bool] = '_vcm',
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
            touched: str | PropertyMap[Vertex, bool] = '_touched',
            faces: Optional[Faces] = None,
            face_patch_map: Optional[PropertyMap, int] = None,
    ) -> None:
        """Perform isotropic remeshing on the specified faces (default: all faces).

        Vertices and edges can be constrained by setting their corresponding values in the
        vertex_constrained and edge_constrained maps to True. Constrained vertices cannot be
        modified during remeshing. Constrained edges *can* be split or collapsed, but not flipped,
        nor its endpoints moved. (If protect_constraints=True, constrained edges cannot be split or
        collapsed.

        If an optional `touched_map: PropertyMap[Vertex, bool]` mapping vertices to bools is
        specified, vertices that were created or moved during the remeshing are flagged as True.
        """
        faces = self.faces if faces is None else faces
        with (
            self.vertex_data.get_or_temp(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm,
            self.vertex_data.get_or_temp(touched, temp_name='_touched', default=False) as touched,
        ):
            args = [
                self.mesh, faces.indices, target_edge_length, n_iter, protect_constraints, vcm.pmap, ecm.pmap, touched.pmap]

            if face_patch_map:
                sgm.meshing.uniform_isotropic_remeshing2(*args, face_patch_map.pmap)
            else:
                sgm.meshing.uniform_isotropic_remeshing(*args)

    def remesh_adaptive(
            self,
            edge_len_min_max: Tuple[float, float],
            tolerance: float,
            n_iter: int = 1,
            ball_radius: float = -1.0,
            protect_constraints=False,
            vertex_constrained: str | PropertyMap[Vertex, bool] = '_vcm',
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
            touched: str | PropertyMap[Vertex, bool] = '_touched',
            faces: Optional[Faces] = None,
            face_patch_map: Optional[PropertyMap, int] = None,
    ) -> None:
        """Isotropic remeshing with a sizing field adaptive to the local curvature.

        Smaller tolerance values lead to shorter edges.

        Ball radius is the radius over which to calculate the local curvature. If ball_radius == -1,
        the 1-ring of each vertex is used.

        For other parameter descriptions see `remesh`.
        """
        faces = self.faces if faces is None else faces
        with (
            self.vertex_data.get_or_temp(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm,
            self.vertex_data.get_or_temp(touched, temp_name='_touched', default=False) as touched,
        ):
            args = [self.mesh, faces, tolerance, ball_radius, edge_len_min_max,
               n_iter, protect_constraints, vcm.pmap, ecm.pmap, touched.pmap]
            if face_patch_map:
                sgm.meshing.adaptive_isotropic_remeshing2(*args, face_patch_map.pmap)
            else:
                sgm.meshing.adaptive_isotropic_remeshing(*args)

    def fair(self, verts: Vertices, continuity=0) -> None:
        """Fair the specified mesh vertices"""
        sgm.meshing.fair(self.mesh, verts, continuity)

    def refine(self, faces: Faces, density=np.sqrt(3)) -> Tuple[Vertices, Faces]:
        """Refine the specified mesh faces

        The number of faces is increased by a factor of `density`.
        Returns indices to the newly created vertices and faces.
        """
        return sgm.meshing.refine(self.mesh, faces, density)

    def smooth_angle_and_area(
            self,
            faces: Faces,
            n_iter: int,
            use_area_smoothing=True,
            use_angle_smoothing=True,
            use_safety_constraints=False,
            do_project=True,
            vertex_constrained: str | PropertyMap[Vertex, bool] = '_vcm',
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
    ) -> None:
        """Smooths a triangulated region of a polygon mesh

        This function attempts to make the triangle angle and area distributions as uniform as possible by moving
        (non-constrained) vertices.
        """
        with (
            self.vertex_data.get_or_temp(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm,
        ):
            sgm.meshing.smooth_angle_and_area(
                self.mesh, faces, n_iter, use_area_smoothing, use_angle_smoothing,
                use_safety_constraints, do_project, vcm.pmap, ecm.pmap)

    def tangential_relaxation(
            self,
            verts: Vertices,
            n_iter: int,
            relax_constraints=False,
            vertex_constrained: str | PropertyMap[Vertex, bool] = '_vcm',
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
    ) -> None:
        with (
            self.vertex_data.get_or_temp(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm,
        ):
            sgm.meshing.tangential_relaxation(
                self.mesh, verts, n_iter, relax_constraints, vcm.pmap, ecm.pmap)

    def smooth_shape(
            self,
            faces: Faces,
            time: float,
            n_iter: int,
            vertex_constrained: str | PropertyMap[Vertex, bool] = '_vcm',
    ) -> None:
        """Smooth the mesh shape by mean curvature flow

        A larger time step results in faster convergence but details may be distorted to a larger extent compared to
         more iterations with a smaller step. Typical values scale in the interval (1e-6, 1]
        """
        with self.vertex_data.get_or_temp(vertex_constrained, temp_name='_vcm', default=False) as vcm:
            sgm.meshing.smooth_shape(self.mesh, faces, time, n_iter, vcm.pmap)

    def does_self_intersect(self) -> bool:
        """Returns True if the mesh self-intersects"""
        return sgm.meshing.does_self_intersect(self.mesh)

    def do_intersect(self, other: Mesh3) -> bool:
        return sgm.meshing.do_intersect(self.mesh, other.mesh)

    def self_intersections(self) -> Tuple[Faces, Faces]:
        """Returns pairs of intersecting faces"""
        return sgm.meshing.self_intersections(self.mesh)

    def remove_self_intersections(self) -> None:
        return sgm.meshing.remove_self_intersections(self.mesh)

    def shortest_path(
            self,
            src_face: Face,
            src_bc: np.ndarray,
            tgt_face: Face,
            tgt_bc: np.ndarray,
    ):
        """Constructs the shortest path between the source and target locations

        locations are specified as a face and barycentric coordinates
        """
        return sgm.locate.shortest_path(self.mesh, src_face, src_bc, tgt_face, tgt_bc)

    def lscm(self, uv_map: str | PropertyMap[Vertex, Point2], initial_verts: Tuple[Vertex, Vertex] = None) -> None:
        """Performs least-squares conformal mapping

        `initial_verts` are indices into the UV map whose coordinates have been fixed.
        Raises a seagullmesh.ParametrizationError if parametrization fails.

        """
        if isinstance(uv_map, str):
            uv_map = self.vertex_data.get(uv_map, default=Point2(0, 0))
        if initial_verts is not None:
            msg = sgm.parametrize.lscm(self.mesh, uv_map.pmap, *initial_verts)
        else:
            msg = sgm.parametrize.lscm(self.mesh, uv_map.pmap)

        if msg != "Success":
            raise ParametrizationError(msg)

    def arap(self, uv_map: str | PropertyMap[Vertex, Point2]) -> None:
        """Performs as-rigid-as-possible parameterization

        Raises a seagullmesh.ParametrizationError if parametrization fails.
        """
        uv_map = self.vertex_data.get(uv_map, default=Point2(0, 0))
        msg = sgm.parametrize.arap(self.mesh, uv_map.pmap)
        if msg != "Success":
            raise ParametrizationError(msg)

    def label_border_vertices(self, is_border: str | PropertyMap[Vertex, bool]):
        is_border = self.vertex_data.get(is_border, default=False)
        sgm.border.label_border_vertices(self.mesh, is_border.pmap)
        return is_border

    def label_border_edges(self, is_border: str | PropertyMap[Edge, bool]):
        is_border = self.edge_data.get(is_border, default=False)
        sgm.border.label_border_edges(self.mesh, is_border.pmap)
        return is_border

    def extract_boundary_cycles(self) -> Halfedges:
        return Halfedges(self, sgm.border.extract_boundary_cycles(self.mesh))

    def has_boundary(self) -> bool:
        return sgm.border.has_boundary(self.mesh)

    def trace_boundary_from_vertex(self, vertex: Vertex) -> Vertices:
        return Vertices(self, sgm.border.trace_boundary_from_vertex(self.mesh, vertex))

    def merge_duplicated_vertices_in_boundary_cycles(self):
        sgm.border.merge_duplicated_vertices_in_boundary_cycles(self.mesh)

    def stitch_borders(self):
        sgm.border.stitch_borders(self.mesh)

    def remesh_planar_patches(
            self,
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
            # face_patch_map: FaceMap = '_face_map',
            cosine_of_maximum_angle: float = 1.0,
    ) -> Mesh3:
        with self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm:
            # fpm = self.face_data.get_or_create_property(face_patch_map, default=-1)
            # TODO see comments in c++ remesh_planar_patches regarding face_patch_map
            out = sgm.meshing.remesh_planar_patches(
                self.mesh, ecm.pmap, cosine_of_maximum_angle)

        return Mesh3(out)

    def edge_collapse(
            self,
            stop_policy_mode: Literal["face", "edge", "edge_length"],
            stop_policy_thresh: float | int,
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
    ) -> int:
        """Mesh simplification by edge collapse

        See https://doc.cgal.org/latest/Surface_mesh_simplification/index.html for details.

        If `stop_policy_mode` is 'face' or 'edge', 'stop_policy_thresh' is either an int or float.
        If an int, stops after the number of edges/faces drops below that threshold. If a float, the
        threshold indicates a ration in [0, 1], stopping when the number of edges/faces is below
        that ratio of the original number.

        If `stop_policy_mode` is 'edge_length', stop_policy_thresh must be a thresh, indicating
        the minimum edge length.

        Constrained edges can be indicated in a boolean edge property map `edge_constrained`.
        """
        if stop_policy_mode == 'face':
            if isinstance(stop_policy_thresh, float):
                fn = sgm.simplification.edge_collapse_face_count_ratio
            elif isinstance(stop_policy_thresh, int):
                fn = sgm.simplification.edge_collapse_face_count
            else:
                raise ValueError(f"Unsupported threshold type {type(stop_policy_thresh)}")
        elif stop_policy_mode == 'edge':
            if isinstance(stop_policy_thresh, float):
                fn = sgm.simplification.edge_collapse_edge_count_ratio
            elif isinstance(stop_policy_thresh, int):
                fn = sgm.simplification.edge_collapse_edge_count
            else:
                raise ValueError(f"Unsupported threshold type {type(stop_policy_thresh)}")
        elif stop_policy_mode == 'edge_length':
            fn = sgm.simplification.edge_collapse_edge_length
        else:
            raise ValueError(f"Unsupported stop policy mode {stop_policy_mode}")

        with self.edge_data.get_or_temp(edge_constrained, temp_name='_ecm', default=False) as ecm:
            out = fn(self.mesh, stop_policy_thresh, ecm.pmap)

        return out

    def skeletonize(self):
        """Construct the medial axis skeleton

        From [Triangulated Surface Mesh Skeletonization](
            https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html).

        Returns a `Skeleton` object with properties
          points : (n, 3) array of medial axis vertex positions
          edges : (m, 2) array of vertex indices
          vertex_map : dict[int, list[mesh vertex]] mapping skeleton vertices to mesh vertices
        """
        from seagullmesh.skeleton import Skeleton
        skel = sgm.skeletonization.extract_mean_curvature_flow_skeleton(self.mesh)
        return Skeleton(mesh=self, skeleton=skel)

    def interpolated_corrected_curvatures(
            self,
            ball_radius: float = -1,
            mean_curvature_map: str | PropertyMap[Vertex, float] = 'mean_curvature',
            gaussian_curvature_map: str | PropertyMap[Vertex, float] = 'gaussian_curvature',
            principal_curvature_map: str | PropertyMap[Vertex, sgm.properties.PrincipalCurvaturesAndDirections] = 'principal_curvature',
    ):
        mcm = self.vertex_data.get(mean_curvature_map, 0.0)
        gcm = self.vertex_data.get(gaussian_curvature_map, 0.0)
        pcm = self.vertex_data.get(
            principal_curvature_map, sgm.properties.PrincipalCurvaturesAndDirections())

        sgm.meshing.interpolated_corrected_curvatures(
            self.mesh, mcm.pmap, gcm.pmap, pcm.pmap, ball_radius)

    @staticmethod
    def from_poisson_surface_reconstruction(
            points: np.ndarray,
            normals: np.ndarray,
            spacing: float,
    ) -> Mesh3:
        mesh = sgm.poisson_reconstruct.reconstruct_surface(points, normals, spacing)
        return Mesh3(mesh)

    def alpha_wrapping(
            self,
            alpha: float | None,
            offset: float | None,
            relative_alpha: float = 20,
            relative_offset: float = 600,
    ) -> Mesh3:
        if alpha is None or offset is None:
            diagonal = self.bounding_box().diagonal()
            alpha = diagonal / relative_alpha if alpha is None else alpha
            offset = diagonal / relative_offset if offset is None else offset

        mesh = sgm.alpha_wrapping.wrap_mesh(self.mesh, alpha, offset)
        return Mesh3(mesh)

    @staticmethod
    def from_alpha_wrapping(
            points: np.ndarray,
            alpha: float | None,
            offset: float | None,
            relative_alpha: float = 20,
            relative_offset: float = 600,
    ):
        if alpha is None or offset is None:
            diagonal = _bbox_diagonal(points)
            alpha = diagonal / relative_alpha if alpha is None else alpha
            offset = diagonal / relative_offset if offset is None else offset

        mesh = sgm.alpha_wrapping.wrap_points(points, alpha, offset)
        return Mesh3(mesh)

    def triangulate_faces(self, faces: Faces | None = None) -> Mesh3:
        if faces:
            sgm.geometry.triangulate_faces(self.mesh, faces)
        else:
            sgm.geometry.triangulate_faces(self.mesh)
        return self

    def reverse_face_orientation(self, faces: Faces | None = None) -> Mesh3:
        if faces:
            sgm.geometry.reverse_face_orientations(self.mesh, faces)
        else:
            sgm.geometry.reverse_face_orientations(self.mesh)
        return self

    def regularize_face_selection_borders(
            self,
            is_selected: PropertyMap[Face, bool],
            weight: float,
            prevent_unselection: bool = False,
    ):
        sgm.border.regularize_face_selection_borders(self.mesh, is_selected, weight, prevent_unselection)

    def label_selected_face_patches(self, faces: Faces, face_patch_idx: PropertyMap[Face, int] | str):
        # faces not in faces are labeled face_patch_idx=0, otherwise 1 + the index of the patch of selected regions
        face_patch_idx = self.face_data.get(face_patch_idx, default=0, is_index=True)
        sgm.connected.label_selected_face_patches(self.mesh, faces, face_patch_idx.pmap)
        return face_patch_idx

    def label_connected_components(
            self,
            face_patches: str | PropertyMap[Face, int] = '_face_patch',
            edge_is_constrained: str | PropertyMap[Edge, bool] ='_ecm',
    ) -> int:
        # Returns number of components

        with (
            self.face_data.get_or_temp(face_patches, '_face_patch', default=0, dtype='uint32') as face_patches,
            self.edge_data.get_or_temp(edge_is_constrained, '_ecm', default=False) as ecm,
        ):
            return sgm.connected.label_connected_components(self.mesh, face_patches.pmap, ecm.pmap)

    def remove_connected_face_patches(
            self,
            to_remove: Sequence[int | bool],
            face_patches: str | PropertyMap[Face, int | bool],
            inplace=False,
    ):
        if not inplace and not isinstance(face_patches, str):
            raise TypeError("Can't supply a pre-existing property map if the mesh is to be copied")

        out = self if inplace else self.copy()
        face_patches = out.face_data.check(face_patches)
        sgm.connected.remove_connected_face_patches(out.mesh, to_remove, face_patches.pmap)
        return out

    def keep_connected_face_patches(
            self,
            to_keep: Sequence[int | bool],
            face_patches: str | PropertyMap[Face, int | bool],
            inplace=False,
    ):
        out = self if inplace else self.copy()
        face_patches = self.face_data.check(face_patches)
        sgm.connected.keep_connected_face_patches(self.mesh, to_keep, face_patches.pmap)
        return out

    def connected_component(
            self, seed_face: Face, edge_is_constrained: PropertyMap[Edge, bool] | str = '_ecm') -> Faces:
        with self.face_data.get_or_temp(edge_is_constrained, tempname='_ecm', default=False) as ecm:
            return sgm.connected.connected_component(self.mesh, seed_face, ecm.pmap)


class ParametrizationError(RuntimeError):
    pass


def _bbox_diagonal(points: np.ndarray):
    x0, y0, z0 = points.min(axis=0)
    x1, y1, z1 = points.max(axis=0)
    return np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)


Key = TypeVar('Key', Vertex, Face, Edge, Halfedge)
Val = TypeVar('Val', int, bool, float, Point2, Point3, Vector2, Vector3)

"""Something that can vector-index a property map"""
IntoIndices = np.ndarray | slice | Sequence[Key]


class PropertyMap(Generic[Key, Val], ABC):
    def __init__(self, pmap, data: MeshData[Key], dtype: str):
        self.pmap = pmap  # the C++ object
        self.data = data
        self._dtype_name = dtype

    def __str__(self) -> str:
        return f'PropertyMap[{self.data.key_type.__name__}, {self.dtype_name}]'

    @property
    def dtype_name(self) -> str:
        return self._dtype_name

    def _check_cpp_indices(self, idxs: _CppIndicesTypes) -> _CppIndicesTypes:
        if isinstance(idxs, self.data.cpp_indices_type):
            return idxs

        msg = f"Tried to index a {self.data.key_type} property map with {idxs} indices"
        raise TypeError(msg)

    def _to_cpp_indices(self, key: IntoIndices[Key]) -> _CppIndicesTypes:
        # Returns the C++ indexer -- sgm.mesh.Vertices, sgm.mesh.Faces, etc
        if isinstance(key, Indices):
            if key.index_type is self.data.key_type:
                return key.indices
            else:
                msg = f'Tried to index a {self.data.key_type} property map with {key.index_type} indices'
                raise TypeError(msg)
        elif isinstance(key, _CppIndicesTypes):  # This should only happen internally..
            return self._check_cpp_indices(key)

        # Two remaining possibilities:
        try:
            # 1) list[key] or similar
            # Free function std::vector<Index> -> Indices<Index>
            idxs = sgm.mesh.make_indices(key)
            return self._check_cpp_indices(idxs)
        except TypeError:
            # 2) index into indices, like a ndarray[bool]
            all_idxs = self.data.py_indices_type.all_indices(self.data.mesh)
            return all_idxs[key].indices

    @overload
    def __getitem__(self, key: int | Key) -> Val: ...

    @overload
    def __getitem__(self, key: IntoIndices[Key]) -> np.ndarray: ...

    def __getitem__(self, key):
        if isinstance(key, int):  # e.g. pmap[0] -> value for the first face
            return self.pmap[self.data.py_indices_type.all_indices(self.data.mesh)[key]]
        elif isinstance(key, _IndexTypes):  # e.g. pmap[Face] -> scalar value
            if isinstance(key, self.data.key_type):
                return self.pmap[key]
            else:
                msg = f'Tried to index a {self.data.key_type} property map with {type(key)} index'
                raise TypeError(msg)
        else:
            return self.pmap[self._to_cpp_indices(key)]

    @overload
    def __setitem__(self, key: int | Key, val: Val) -> None: ...

    @overload
    def __setitem__(self, key: IntoIndices[Key], val: np.ndarray | Sequence[Val]) -> None: ...

    def __setitem__(self, key, val) -> None:
        if isinstance(key, int):  # pmap[0] -> value for the first face
            self.pmap[self.data.all_mesh_keys[key]] = val
        elif isinstance(key, _IndexTypes):
            if isinstance(key, self.data.key_type):
                self.pmap[key] = val
            else:
                msg = f'Tried to index a {self.data.key_type} property map with {type(key)} index'
                raise TypeError(msg)
        else:
            self.pmap[self._to_cpp_indices(key)] = val

    for dunder in (
            '__add__',
            '__eq__',
            '__ge__',
            '__gt__',
            '__le__',
            '__lt__',
            '__mul__',
            '__ne__',
            '__neg__',
            '__pos__',
            '__pow__',
            '__mod__',
    ):
        def _dunder_impl(self, other, _dunder=dunder):
            if isinstance(other, PropertyMap):
                other = other[:]
            fn = getattr(self[:], _dunder)
            return fn(other)

        locals()[dunder] = _dunder_impl


class ScalarPropertyMap(PropertyMap[Key, Val]):
    pass


class ArrayPropertyMap(PropertyMap[Key, Val]):
    def get_objects(self, key: IntoIndices[Key]) -> Sequence[Val]:
        # __getitem__ defaults to returning array(nk, ndim), also allow returning list[Point2]
        return self.pmap.get_vector(self._to_cpp_indices(key))

    def set_objects(self, key: IntoIndices[Key], val: Sequence[Val]) -> None:
        # handled by the general case
        self[key] = val


_PMapDType = str | np.dtype | type
VertexIndexMap = PropertyMap[Vertex, int]
VertexPointMap = PropertyMap[Vertex, Point2 | Point3]


class MeshData(Generic[Key]):
    py_indices_type: Type[_PyIndicesTypes]  # set by subclass

    # set in __init_subclass__
    key_type: Type[Key]  # == py_indices_type.index_type
    cpp_indices_type: Type[_CppIndicesTypes]
    _key_name: str  # 'vertices', 'faces', 'edges', 'halfedges'
    _prefix: str  # V, F, E, H

    def __init__(self, mesh: Mesh3):
        self._data: Dict[str, PropertyMap[Key]] = {}
        self.mesh = mesh  # python wrapped mesh

    _dtype_mappings: dict[type, str] = {
        float: 'double',
        int: 'int64',
    }

    @property
    def all_mesh_keys(self) -> Indices:
        return self.py_indices_type.all_indices(self.mesh)

    @property
    def n_mesh_keys(self) -> int:
        return self.py_indices_type.n_mesh_keys(self.mesh)

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.key_type = cls.py_indices_type.index_type
        cls.cpp_indices_type = cls.py_indices_type.indices_type
        cls._key_name = cls.py_indices_type.__name__.lower()  # 'vertices', 'faces', etc
        cls._prefix = cls._key_name[0].upper()  # 'V', 'F', 'E', 'H'

    def __repr__(self) -> str:
        pmaps = ', '.join(f'{k}: {v}' for k, v in self.items())
        return f'{self.__class__.__name__}[{self.key_type.__name__}]({self.mesh}): {pmaps}'

    def _dtype_name(self, dtype: _PMapDType) -> str:
        if isinstance(dtype, str):
            return dtype
        elif isinstance(dtype, np.dtype):
            if dtype.name == 'float64':
                return 'double'
            return dtype.name
        elif mapped := self._dtype_mappings.get(dtype):
            return mapped
        else:
            return dtype.__name__

    def _pmap_class_name(self, dtype_name: str) -> str:
        return f'{self._prefix}_{dtype_name}_PropertyMap'

    @contextmanager
    def temp(
            self,
            name: str,
            default: Val,
            dtype: _PMapDType | None = None,
    ) -> Iterator[PropertyMap[Key, Val]]:
        """Create a temporary property map and remove it after exiting the contextmanager"""
        pmap = self.create(name=name, default=default, dtype=dtype)
        yield pmap
        self.remove(name)

    @contextmanager
    def get_or_temp(
            self,
            pmap: str | PropertyMap[Key, Val] | None,
            temp_name: str,
            default: Val,
            dtype: _PMapDType | None = None,
    ) -> Iterator[PropertyMap[Key, Val]]:
        """Construct a placeholder property map

        If pmap: str == temp_name, the map is removed after use. If pmap is another str or a
        pre-existing property map, it's not removed.
        """
        remove_after = isinstance(pmap, str) and pmap == temp_name
        pmap = self.get(pmap or temp_name, default=default, dtype=dtype)
        yield pmap
        if remove_after:
            self.remove(temp_name)

    def create(
            self,
            name: str,
            default: Val,
            dtype: str | np.dtype | None = None,
    ) -> PropertyMap[Key, Val]:
        """Add a property map

        The type of the map's value is inferred from the default value. E.g.
        `mesh.vertex_data.add_property('foo', 0.0)` constructs a property map named 'foo'
        storing doubles on vertices and `mesh.face_data.add_property('bar', Point2(0, 0))` stores
        2d points on mesh faces.

        Mapping int-values to C++ types is inherently ambiguous, so it's preferred to explicitly
        specify a dtype, either as a string like 'uint32' or a numpy dtype.
        """
        dtype_name = self._dtype_name(dtype or type(default))
        cls_name = self._pmap_class_name(dtype_name)

        try:
            pmap_class: type = getattr(sgm.properties, cls_name)
        except AttributeError:
            msg = (
                f"Property map class {cls_name} does not exist for default = "
                f"{type(default)}({default}) with supplied {dtype=} and inferred {dtype_name=} "
                f"during attempted construction of property name {name=}"
            )
            raise TypeError(msg)

        cpp_pmap = pmap_class(self.mesh.mesh, name, default)
        wrapped_pmap = self._data[name] = self.wrap(cpp_pmap=cpp_pmap, dtype_name=dtype_name)
        return wrapped_pmap

    def remove(self, key: str):
        pmap = self._data.pop(key)
        sgm.properties.remove_property_map(self.mesh.mesh, pmap.pmap)

    def find(
            self,
            pmap_cls: type,
            name: str,
            wrapper_cls: Type[PropertyMap] | None = None,
            dtype_name: str = 'unknown',
    ) -> PropertyMap[Key]:
        # Finds a pre-existing pmap, wraps it, and returns it without adding it to self
        pmap = pmap_cls.get_property_map(self.mesh.mesh, name)  # noqa
        if pmap is None:
            raise KeyError(f"Property map {pmap_cls} {name} doesn't exist")
        return self.wrap(pmap, wrapper_cls=wrapper_cls, dtype_name=dtype_name)

    def check(self, item: str | PropertyMap[Key, Val]) -> PropertyMap[Key, Val]:
        pmap = self[item] if isinstance(item, str) else item
        return self._check_property_map(pmap)

    def wrap(
            self,
            cpp_pmap,  # the c++ class,
            wrapper_cls: Type[PropertyMap] | None = None,
            dtype_name: str = 'unknown',
    ) -> PropertyMap[Key]:
        if wrapper_cls is None:
            wrapper_cls = ScalarPropertyMap if cpp_pmap.is_scalar else ArrayPropertyMap
        return wrapper_cls(pmap=cpp_pmap, data=self, dtype=dtype_name)

    def _check_property_map(self, pmap: PropertyMap) -> PropertyMap[Key, Any]:
        if pmap.data is not self:
            msg = f'Trying to get {pmap} with mesh {pmap.data.mesh} from {self} with mesh {self.mesh}'
            raise TypeError(msg)

        return pmap

    def get(
            self,
            key: str | PropertyMap[Key, Val],
            default: Val,
            dtype: _PMapDType | None = None,
    ) -> PropertyMap[Key, Val]:
        if isinstance(key, PropertyMap):
            return self._check_property_map(key)

        if key in self._data:
            return self._data[key]
        else:
            assert isinstance(key, str), key
            return self.create(name=key, default=default, dtype=dtype)

    def __getitem__(self, item: str) -> PropertyMap[Key, Any]:
        return self._data[item]

    def __delitem__(self, item: str):
        self.remove(item)

    def __setitem__(self, key: str, value: Any) -> None:
        if isinstance(value, PropertyMap):  # An already wrapped pmap
            if value.data is not self:
                raise TypeError("Property map does not belong to this mesh")
            self._data[key] = value
            return

        if hasattr(value, "is_sgm_property_map"):
            # Assigning the bare C++ property map
            self._data[key] = self.wrap(cpp_pmap=value)
            return

        # Implicit construction of a new property map with initial value(s) `value`
        default = np.zeros_like(value, shape=()).item()  # type: ignore
        pmap = self.get(key, default)
        pmap[self.py_indices_type.all_indices(self.mesh)] = value

    def items(self) -> Iterator[Tuple[str, PropertyMap[Key, Any]]]:
        yield from self._data.items()

    def values(self) -> Iterator[PropertyMap[Key, Any]]:
        yield from self._data.values()

    def keys(self) -> Iterator[str]:
        yield from self._data.keys()

    def __iter__(self) -> Iterator[str]:
        yield from self._data.__iter__()

    def check_has_same_properties(self, other: Self) -> None:
        if missing_keys := set(self.keys()).symmetric_difference(other.keys()):
            msg = f"{self.key_type} properties {missing_keys} are not present in both meshes."  # noqa
            raise ValueError(msg)

        for k, pmap in self.items():
            t0, t1 = type(pmap.pmap), type(other[k].pmap)
            if t0 is not t1:
                raise TypeError(f"Property {k} has two different types: {t0} and {t1}")


class VertexData(MeshData):
    py_indices_type = Vertices


class FaceData(MeshData):
    py_indices_type = Faces


class EdgeData(MeshData):
    py_indices_type = Edges


class HalfedgeData(MeshData):
    py_indices_type = Halfedges
