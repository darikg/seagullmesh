from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from contextlib import contextmanager
from functools import cached_property
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING, Union, Sequence, TypeVar, overload, Tuple, \
    Generic, List, Iterator, Type, Dict, Literal

import numpy as np

from seagullmesh._seagullmesh.mesh import (
    Mesh3 as _Mesh3,
    Point2, Point3, Vector2, Vector3,
    Vertex, Face, Edge, Halfedge,
)
from seagullmesh import _seagullmesh as sgm
from ._version import version_info, __version__  # noqa

Vertices = Sequence[Vertex]
Faces = Sequence[Face]
Edges = Sequence[Edge]
Halfedges = Sequence[Halfedge]

if TYPE_CHECKING:
    try:
        import pyvista as pv  # noqa
    except ImportError:
        pv = None


class ParametrizationError(RuntimeError):
    pass


class Mesh3:
    def __init__(self, mesh: _Mesh3 | None = None):
        if mesh is None:
            mesh = _Mesh3()

        self._mesh = mesh

        if hasattr(sgm, 'properties'):
            self.vertex_data = MeshData(mesh, sgm.properties.add_vertex_property, 'vertices', 'Vert')
            self.face_data = MeshData(mesh, sgm.properties.add_face_property, 'faces', 'Face')
            self.edge_data = MeshData(mesh, sgm.properties.add_edge_property, 'edges', 'Edge')
            self.halfedge_data = MeshData(mesh, sgm.properties.add_halfedge_property, 'halfedges', 'HalfEdge')
        else:
            warnings.warn("properties module not available")

        # Mesh3 automatically constructs a vertex point 3 property, make it available
        # from the python side
        # self.vertex_data.assign_property_map('points', mesh.points)

    @property
    def vertices(self) -> sgm.mesh.Vertices:
        """Vector of vertex indices"""
        return self._mesh.vertices

    @property
    def faces(self) -> sgm.mesh.Faces:
        """Vector of face indices"""
        return self._mesh.faces

    @property
    def edges(self) -> sgm.mesh.Edges:
        """Vector of edge indices"""
        return self._mesh.edges

    @property
    def halfedges(self) -> sgm.mesh.Halfedges:
        """Vector of halfedge indices"""
        return self._mesh.halfedges

    n_vertices = property(lambda self: self._mesh.n_vertices)
    n_faces = property(lambda self: self._mesh.n_faces)
    n_edges = property(lambda self: self._mesh.n_edges)
    n_halfedges = property(lambda self: self._mesh.n_halfedges)

    is_valid = property(lambda self: self._mesh.is_valid)

    @property
    def mesh(self) -> sgm.mesh.Mesh3:
        """The C++ mesh object"""
        return self._mesh

    def copy(self) -> Mesh3:
        """Deep-copy the mesh and all its properties"""
        out = Mesh3(sgm.mesh.Mesh3(self._mesh))

        # The properties have been copied, just need to create new pmap references
        for k in ('vertex_data', 'face_data', 'edge_data', 'halfedge_data'):
            _copy_property_metadata(getattr(self, k), out._mesh, getattr(out, k))

        return out

    def transform(self, transform: np.ndarray, inplace=False) -> Mesh3:
        out = self if inplace else self.copy()
        out._mesh.transform(transform)
        return out

    def add(self, other: Mesh3, check_properties=True):
        pass

    @property
    def has_garbage(self) -> bool:
        return self._mesh.has_garbage

    def collect_garbage(self) -> None:
        self._mesh.collect_garbage()

    def volume(self) -> float:
        return self._mesh.volume()

    def expand_selection(self, selection: Sequence[Key]) -> Sequence[Key]:
        """Given a list of vertices or faces, returns a sequence containing the original and adjacent elements"""
        return self._mesh.expand_selection(selection)

    def vertices_to_faces(self, verts: Vertices) -> Faces:
        return self._mesh.vertices_to_faces(verts)

    def to_polygon_soup(self) -> Tuple[np.ndarray, np.ndarray]:
        """Returns vertices (nv * 3) and faces (nf * 3) array"""
        return self._mesh.to_polygon_soup()

    def face_normals(self, faces: Faces) -> np.ndarray:
        """Returns a (len(faces) * 3) array of face normal vectors"""
        return self._mesh.face_normals(faces)

    def bounding_box(self) -> sgm.mesh.BoundingBox3:
        return self._mesh.bounding_box()

    @staticmethod
    def from_polygon_soup(verts: np.ndarray, faces: np.ndarray, orient=True) -> Mesh3:
        """Constructs a surface mesh from vertices (nv * 3) and faces (nf * 3) arrays

        If `orient` is True (default), the faces are reindexed to represent a consistent manifold surface.
        """
        mesh = sgm.mesh.polygon_soup_to_mesh3(verts, faces, orient)
        return Mesh3(mesh)

    @staticmethod
    def from_file(filename: str) -> Mesh3:
        mesh = sgm.mesh.load_mesh_from_file(filename)
        return Mesh3(mesh)

    def to_file(self, filename: str):
        ext = Path(filename).suffix
        if ext == '.ply':
            self._mesh.write_ply(filename)
        elif ext == '.off':
            self._mesh.write_off(filename)
        else:
            raise ValueError(f"Unsupported format '{ext}'")

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
            vertex_data: Literal['all'] | Sequence[str] = (),
            face_data: Literal['all'] | Sequence[str] = (),
    ) -> pv.PolyData:
        """Returns the mesh as a `pyvista.PolyData` object.

        By default, vertex and cell data is ignored -- specify vertex_data and cell_data as a list of keys
        naming property maps to copy, or the string 'all' for all of them.
        """
        import pyvista as pv
        verts, faces = self._mesh.to_polygon_soup()
        mesh = pv.PolyData.from_regular_faces(verts, faces)

        if vertex_data:
            keys = self.vertex_data.keys() if vertex_data == 'all' else vertex_data
            vertices = self.vertices
            for k in keys:
                mesh.point_data[k] = self.vertex_data[k][vertices]

        if face_data:
            keys = self.face_data.keys() if face_data == 'all' else face_data
            faces = self.faces
            for k in keys:
                mesh.cell_data[k] = self.face_data[k][faces]

        return mesh

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
            self.vertex_data.temporary(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm,
            self.vertex_data.temporary(touched, temp_name='_touched', default=False) as touched,
        ):
            args = [
                self._mesh, faces, target_edge_length, n_iter, protect_constraints, vcm.pmap, ecm.pmap, touched.pmap]

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
            self.vertex_data.temporary(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm,
            self.vertex_data.temporary(touched, temp_name='_touched', default=False) as touched,
        ):
            args = [self._mesh, faces, tolerance, ball_radius, edge_len_min_max,
               n_iter, protect_constraints, vcm.pmap, ecm.pmap, touched.pmap]
            if face_patch_map:
                sgm.meshing.adaptive_isotropic_remeshing2(*args, face_patch_map.pmap)
            else:
                sgm.meshing.adaptive_isotropic_remeshing(*args)

    def fair(self, verts: Vertices, continuity=0) -> None:
        """Fair the specified mesh vertices"""
        sgm.meshing.fair(self._mesh, verts, continuity)

    def refine(self, faces: Faces, density=np.sqrt(3)) -> Tuple[Vertices, Faces]:
        """Refine the specified mesh faces

        The number of faces is increased by a factor of `density`.
        Returns indices to the newly created vertices and faces.
        """
        return sgm.meshing.refine(self._mesh, faces, density)

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
            self.vertex_data.temporary(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm,
        ):
            sgm.meshing.smooth_angle_and_area(
                self._mesh, faces, n_iter, use_area_smoothing, use_angle_smoothing,
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
            self.vertex_data.temporary(vertex_constrained, temp_name='_vcm', default=False) as vcm,
            self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm,
        ):
            sgm.meshing.tangential_relaxation(
                self._mesh, verts, n_iter, relax_constraints, vcm.pmap, ecm.pmap)

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
        with self.vertex_data.temporary(vertex_constrained, temp_name='_vcm', default=False) as vcm:
            sgm.meshing.smooth_shape(self._mesh, faces, time, n_iter, vcm.pmap)

    def does_self_intersect(self) -> bool:
        """Returns True if the mesh self-intersects"""
        return sgm.meshing.does_self_intersect(self._mesh)

    def self_intersections(self) -> Tuple[Faces, Faces]:
        """Returns pairs of intersecting faces"""
        return sgm.meshing.self_intersections(self._mesh)

    def remove_self_intersections(self) -> None:
        return sgm.meshing.remove_self_intersections(self._mesh)

    def _get_vertex_point_map(self, vert_points: str | PropertyMap[Vertex, Point2 | Point3] | None = None):
        if vert_points is None:
            return self.mesh.points
        else:
            return self.vertex_data.get_property_map(vert_points).pmap

    def aabb_tree(self, vert_points: str | PropertyMap[Vertex, Point2 | Point3] | None = None):
        """Construct an axis-aligned bounding box tree for accelerated point location by `Mesh3.locate_points

        By default, the AABB tree is constructed for the default mesh vertex locations, but also accepts a vertex
        property map storing Point2 or Point3 locations.
        """
        return sgm.locate.aabb_tree(self._mesh, self._get_vertex_point_map(vert_points))

    def locate_points(
            self,
            points: np.ndarray,
            aabb_tree=None,
            vert_points: str | PropertyMap[Vertex, Point2 | Point3] | None = None,
    ) -> Tuple[Faces, np.ndarray]:
        """Given an array of points, locate the nearest corresponding points on the mesh

        `aabb_tree` is an optional axis-aligned bounding box from Mesh3.aabb_tree. If the tree was constructed with a
        Point2 vertex property map, `points` is of shape (np, 2), otherwise (np, 3).

        Returns a list of face indices of length np, and an array (np, 3) of barycentric coordinates within those faces.
        """
        tree = aabb_tree or self.aabb_tree()
        pmap = self._get_vertex_point_map(vert_points)
        surface_points = sgm.locate.locate_points(self._mesh, tree, points, pmap)
        return surface_points.faces, surface_points.bary_coords

    def construct_points(
            self,
            faces: Faces,
            bary_coords: np.ndarray,
            vert_points: str | PropertyMap[sgm.mesh.Vertex, sgm.mesh.Point2 | sgm.mesh.Point3] | None = None,
    ) -> np.ndarray:
        """Construct a set of points from face barycentric coordinates

        `bary_coords` must be of shape (len(faces), 3)

        By default, points are constructed using the default mesh vertex points. An optional vertex point map of value
        Point2 or Point3 can also be supplied. The returned array if of shape (len(faces), 2 or 3) as appropriate.
        """
        pmap = self._get_vertex_point_map(vert_points)
        return sgm.locate.construct_points(self._mesh, faces, bary_coords, pmap)

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
        return sgm.locate.shortest_path(self._mesh, src_face, src_bc, tgt_face, tgt_bc)

    def lscm(self, uv_map: str | PropertyMap[Vertex, Point2], initial_verts: Tuple[Vertex, Vertex] = None) -> None:
        """Performs least-squares conformal mapping

        `initial_verts` are indices into the UV map whose coordinates have been fixed.
        Raises a seagullmesh.ParametrizationError if parametrization fails.

        """
        if isinstance(uv_map, str):
            uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        if initial_verts is not None:
            msg = sgm.parametrize.lscm(self._mesh, uv_map.pmap, *initial_verts)
        else:
            msg = sgm.parametrize.lscm(self._mesh, uv_map.pmap)

        if msg != "Success":
            raise ParametrizationError(msg)

    def arap(self, uv_map: str | PropertyMap[Vertex, Point2]) -> None:
        """Performs as-rigid-as-possible parameterization

        Raises a seagullmesh.ParametrizationError if parametrization fails.
        """
        uv_map = self.vertex_data.get_or_create_property(uv_map, default=Point2(0, 0))
        msg = sgm.parametrize.arap(self._mesh, uv_map.pmap)
        if msg != "Success":
            raise ParametrizationError(msg)

    def estimate_geodesic_distances(
            self,
            src: Union[Vertex, Vertices],
            distance_prop: str | PropertyMap[Vertex, float],
    ):
        """Estimates the geodesic distance from the source vertex/vertices to all vertices in the mesh

        Estimated distances are stored in the supplied vertex property map.
        """
        distances = self.vertex_data.get_or_create_property(distance_prop, default=0.0)
        self._mesh.estimate_geodesic_distances(distances.pmap, src)

    def label_border_vertices(self, is_border: str | PropertyMap[Vertex, bool]):
        is_border = self.vertex_data.get_or_create_property(is_border, default=False)
        sgm.border.label_border_vertices(self._mesh, is_border.pmap)
        return is_border

    def label_border_edges(self, is_border: str | PropertyMap[Edge, bool]):
        is_border = self.edge_data.get_or_create_property(is_border, default=False)
        sgm.border.label_border_edges(self._mesh, is_border.pmap)
        return is_border

    def extract_boundary_cycles(self) -> Halfedges:
        return sgm.border.extract_boundary_cycles(self._mesh)

    def remesh_planar_patches(
            self,
            edge_constrained: str | PropertyMap[Edge, bool] = '_ecm',
            # face_patch_map: FaceMap = '_face_map',
            cosine_of_maximum_angle: float = 1.0,
    ) -> Mesh3:
        with self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm:
            # fpm = self.face_data.get_or_create_property(face_patch_map, default=-1)
            # TODO see comments in c++ remesh_planar_patches regarding face_patch_map
            out = sgm.meshing.remesh_planar_patches(
                self._mesh, ecm.pmap, cosine_of_maximum_angle)

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

        with self.edge_data.temporary(edge_constrained, temp_name='_ecm', default=False) as ecm:
            out = fn(self._mesh, stop_policy_thresh, ecm.pmap)

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
        skeleton = sgm.skeletonization.extract_mean_curvature_flow_skeleton(self._mesh)
        return Skeleton(mesh=self, skeleton=skeleton)

    def interpolated_corrected_curvatures(
            self,
            ball_radius: float = -1,
            mean_curvature_map: str | PropertyMap[Vertex, float] = 'mean_curvature',
            gaussian_curvature_map: str | PropertyMap[Vertex, float] = 'gaussian_curvature',
            principal_curvature_map: str | VertexPrincipalCurvaturesMap = 'principal_curvature',
    ):
        mcm = self.vertex_data.get_or_create_property(mean_curvature_map, 0.0)
        gcm = self.vertex_data.get_or_create_property(gaussian_curvature_map, 0.0)
        pcm = self.vertex_data.get_or_create_property(
            principal_curvature_map, sgm.properties.PrincipalCurvaturesAndDirections())

        sgm.meshing.interpolated_corrected_curvatures(
            self._mesh, mcm.pmap, gcm.pmap, pcm.pmap, ball_radius)

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

        mesh = sgm.alpha_wrapping.wrap_mesh(self._mesh, alpha, offset)
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

    def triangulate_faces(self, faces: Faces | None = None):
        faces = self.faces if faces is None else faces
        sgm.triangulate.triangulate_faces(self._mesh, faces)

    def reverse_face_orientation(self, faces: Faces | None = None):
        faces = self.faces if faces is None else faces
        sgm.triangulate.reverse_face_orientations(self._mesh, faces)

    def does_bound_a_volume(self) -> bool:
        return sgm.triangulate.does_bound_a_volume(self._mesh)

    def is_outward_oriented(self) -> bool:
        return sgm.triangulate.is_outward_oriented(self._mesh)

    def regularize_face_selection_borders(
            self,
            is_selected: PropertyMap[Face, bool],
            weight: float,
            prevent_unselection: bool = False,
    ):
        sgm.border.regularize_face_selection_borders(self._mesh, is_selected, weight, prevent_unselection)

    def vertex_degrees(self, vertices: Vertices | None = None) -> np.ndarray:
        vertices = self.vertices if vertices is None else vertices
        return sgm.mesh.vertex_degrees(vertices)

    def label_selected_face_patches(self, faces: Faces, face_patch_idx: PropertyMap[Face, int] | str):
        # faces not in faces are labeled face_patch_idx=0, otherwise 1 + the index of the patch of selected regions
        face_patch_idx = self.face_data.get_or_create_property(face_patch_idx, default=0, is_index=True)
        sgm.connected.label_selected_face_patches(self._mesh, faces, face_patch_idx.pmap)
        return face_patch_idx

    def label_connected_components(self, face_patches: PropertyMap[Face, int], edge_is_constrained: PropertyMap[Edge, bool]) -> int:
        return sgm.connected.label_connected_components(self._mesh, face_patches.pmap, edge_is_constrained.pmap)

    def remove_connected_face_patches(self, to_remove: Sequence[int], face_patches: PropertyMap[Face, int]):
        sgm.connected.remove_connected_face_patches(self._mesh, to_remove, face_patches.pmap)

    def connected_component(
            self, seed_face: Face, edge_is_constrained: PropertyMap[Edge, bool] | str = '_ecm') -> Faces:
        with self.face_data.temporary(edge_is_constrained, tempname='_ecm', default=False) as ecm:
            return sgm.connected.connected_component(self._mesh, seed_face, ecm.pmap)

    def edge_soup(self) -> ndarray:
        return self._mesh.edge_soup()

    def triangle_soup(self) -> ndarray:
        return self._mesh.triangle_soup()


def _bbox_diagonal(points: ndarray):
    x0, y0, z0 = points.min(axis=0)
    x1, y1, z1 = points.max(axis=0)
    return sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)


class Skeleton:
    """Wrapper around the C++ sgm.skeletonization.Skeleton class

    (Which is itself a boost adjacency_list)
    """
    def __init__(self, mesh: Mesh3, skeleton):
        self._mesh = mesh
        self._skeleton = skeleton

    @cached_property
    def points(self) -> ndarray:
        return self._skeleton.points

    @cached_property
    def edges(self) -> ndarray:
        return self._skeleton.edges

    @cached_property
    def vertex_map(self) -> Dict[int, List[Vertex]]:
        return self._skeleton.vertex_map

    @cached_property
    def radii(self) -> ndarray:
        return self._skeleton.compute_radii(self._mesh.mesh)

    def to_pyvista(self):
        import pyvista as pv
        sk_mesh = pv.PolyData()
        sk_mesh.points = self.points
        sk_mesh.lines = pv.CellArray.from_regular_cells(self.edges)
        sk_mesh.point_data['min_radius'] = self.radii[:, 0]
        sk_mesh.point_data['max_radius'] = self.radii[:, 1]


Key = TypeVar('Key', Vertex, Face, Edge, Halfedge)
Val = TypeVar('Val', int, bool, float, Point2, Point3, Vector2, Vector3)


class PropertyMap(Generic[Key, Val], ABC):
    def __init__(self, pmap, data: MeshData[Key]):
        self.pmap = pmap  # the C++ object
        self._data = data

    @abstractmethod
    def all_values(self): ...

    @abstractmethod
    def __getitem__(self, key): ...

    @abstractmethod
    def __setitem__(self, key, val): ...

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
                other = other.all_values()
            fn = getattr(self.all_values(), _dunder)
            return fn(other)

        locals()[dunder] = _dunder_impl


class ScalarPropertyMap(PropertyMap[Key, Val]):
    @overload
    def __getitem__(self, key: Union[int, Key]) -> Val: ...

    @overload
    def __getitem__(self, key: Union[A, Sequence[Key], slice]) -> Sequence[Val]: ...

    def __getitem__(self, key):
        try:
            return self.pmap[key]
        except TypeError:
            # Could be some sort of subscripting indexing vector
            return self.pmap[self._data.mesh_keys[key]]

    def __setitem__(self, key, val):
        try:
            self.pmap[key] = val
        except TypeError:
            self.pmap[self._data.mesh_keys[key]] = val

    def all_values(self):
        return self.pmap[self._data.mesh_keys]


class ArrayPropertyMap(PropertyMap[Key, Val]):
    def __getitem__(self, key) -> A:
        try:
            return self.pmap.get_array(key)
        except TypeError:
            return self.pmap.get_array(self._data.mesh_keys[key])

    def __setitem__(self, key, val: A):
        try:
            self.pmap.set_array(key, val)
        except TypeError:
            self.pmap.set_array(self._data.mesh_keys[key], val)

    def get_objects(self, key) -> Sequence[Val]:
        try:
            return self.pmap[key]
        except TypeError:
            return self.pmap[self._data.mesh_keys[key]]

    def set_objects(self, key, val):
        try:
            self.pmap[key] = val
        except TypeError:
            self.pmap[self._data.mesh_keys[key]] = val

    def all_values(self):
        return self.pmap.get_array(self._data.mesh_keys)

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
            if isinstance(other, ScalarPropertyMap):
                other = other.all_values()
            fn = getattr(self.all_values(), _dunder)
            return fn(other)

        locals()[dunder] = _dunder_impl



class MeshData(Generic[Key]):
    def __init__(self, mesh: sgm.mesh.Mesh3, add_fn, key_name: str, prefix: str):
        self._data: Dict[str, PropertyMap[Key]] = {}
        self._mesh = mesh  # the c++ mesh
        self._add_fn = add_fn  # e.g. sgm.properties.add_vertex_property, etc
        self._key_name = key_name  # e.g. 'vertices', 'faces', etc
        self._cpp_cls_prefix = prefix  # e.g. 'Vert', 'Face', etc

    @property
    def mesh_keys(self) -> List[Key]:
        return getattr(self._mesh, self._key_name)

    @property
    def n_mesh_keys(self) -> int:
        return getattr(self._mesh, f'n_{self._key_name}')

    @contextmanager
    def temporary(
            self,
            pmap: str | PropertyMap[Key, Val] | None,
            temp_name: str,
            default: Val,
            signed: bool | None = None,
    ) -> Iterator[PropertyMap[Key, Val]]:
        """Construct a placeholder property map

        If pmap: str == temp_name, the map is removed after use. If pmap is another str or a
        pre-existing property map, it's not removed.
        """
        remove_after = isinstance(pmap, str) and pmap == temp_name
        pmap = self.get_or_create_property(pmap, default=default, signed=signed)
        yield pmap
        if remove_after:
            self.remove_property(temp_name)

    def add_property(
            self,
            key: str,
            default: Val,
            signed: Optional[bool] = None,
            is_index: Optional[bool] = None,
    ) -> PropertyMap[Key, Val]:
        """Add a property map

        The type of the map's value is inferred from the default value. E.g.
        `mesh.vertex_data.add_property('foo', 0.0)` constructs a VertDoublePropertyMap named 'foo'
        storing doubles on vertices and `mesh.face_data.add_property('bar', Point2(0, 0))` stores
        2d points on mesh faces.

        There is an ambiguity when constructing int-valued property maps where some CGAL routines
        require either signed or unsigned ints, so the optional `signed` parameter can make it
        explicit whether signed or unsigned ints are required.
        """
        if signed is is_index is None:
            # Infer C++ property map type from the default
            pmap = self._add_fn(self._mesh, key, default)
        elif isinstance(default, int):
            if is_index:
                # Specify e.g. "FaceIndexPropertyMap"
                pmap_cls = f"{self._cpp_cls_prefix}IndexPropertyMap"
            else:
                assert signed is not None
                # Specify e.g. "VertIntPropertyMap", "VertUIntPropertyMap"
                pmap_cls = f"{self._cpp_cls_prefix}{'' if signed else 'U'}IntPropertyMap"
            pmap = getattr(sgm.properties, pmap_cls)(self._mesh, key, default)
        else:
            raise TypeError(
                f"Can only specify signed for int values, got default={type(default)}({default})")

        return self.assign_property_map(name=key, pmap=pmap)  # The wrapped map

    def remove_property(self, key: str):
        pmap = self._data.pop(key)
        sgm.properties.remove_property_map(self._mesh, pmap.pmap)

    def assign_property_map(
            self,
            name: str,
            pmap,  # The C++ property map
            wrapper_cls: Type[PropertyMap] | None = None
    ) -> PropertyMap:
        if wrapper_cls is None:
            wrapper_cls = ScalarPropertyMap if pmap._is_scalar else ArrayPropertyMap
        wrapped_pmap = self._data[name] = wrapper_cls(pmap=pmap, data=self)
        return wrapped_pmap

    def get_property_map(self, key: str | PropertyMap[Key, Val]) -> PropertyMap[Key, Val]:
        if isinstance(key, PropertyMap):
            return key

        return self._data[key]

    def get_or_create_property(
            self,
            key: str | PropertyMap[Key, Val],
            default: Val,
            signed: bool | None = None,
            is_index: bool | None = None,
    ) -> PropertyMap[Key, Val]:
        if isinstance(key, PropertyMap):
            return key

        if key in self._data:
            return self._data[key]
        else:
            return self.add_property(key, default, signed=signed, is_index=is_index)

    def __getitem__(self, item: str) -> PropertyMap[Key, Any]:
        return self._data[item]

    def __delitem__(self, item: str):
        self.remove_property(item)

    def __setitem__(self, key: str, value: Any):
        if hasattr(value, "_is_sgm_property_map"):
            # Assigning the bare C++ property map
            self.assign_property_map(name=key, pmap=value)
            return

        # Implicit construction of a new property map with initial value(s) `value`
        default = np.zeros_like(value, shape=()).item()
        pmap = self.get_or_create_property(key, default)
        pmap[self.mesh_keys] = value

    def items(self) -> Iterator[Tuple[str, PropertyMap[Key, Any]]]:
        yield from self._data.items()

    def values(self) -> Iterator[PropertyMap[Key, Any]]:
        yield from self._data.values()

    def keys(self) -> Iterator[str]:
        yield from self._data.keys()

    def __iter__(self) -> Iterator[str]:
        yield from self._data.__iter__()


def _copy_property_metadata(src_data: MeshData, dest_mesh: _Mesh3, dest_data: MeshData):
    for name, pmap_wrapper in src_data.items():
        dest_pmap = type(pmap_wrapper.pmap)(dest_mesh, name)
        if dest_pmap is None:
            raise KeyError(f"Property map {name} doesn't exist in destination {type(dest_data)}")
        dest_data.assign_property_map(name, dest_pmap, wrapper_cls=pmap_wrapper.__class__)


class TubeMesher:
    def __init__(self, t0: float, theta0: np.ndarray, pts0: np.ndarray, closed=False):
        mesh = self.mesh = Mesh3()
        t_map = mesh.vertex_data.add_property('t', default=-1.0)
        theta_map = mesh.vertex_data.add_property('theta', default=-1.0)
        is_cap_map = mesh.face_data.add_property('is_cap', default=False)
        self.tube_mesher = sgm.triangulate.TubeMesher(mesh.mesh, t_map.pmap, theta_map.pmap, is_cap_map.pmap, t0, theta0, pts0)

        self.closed = closed
        if closed:
            self.tube_mesher.close_xs(False)

    def add_xs(self, t: float, theta: np.ndarray, pts: np.ndarray):
        self.tube_mesher.add_xs(t, theta, pts)

    def finish(self, reverse_orientation: bool):
        if self.closed:
            self.tube_mesher.close_xs(True)

        sgm.triangulate.triangulate_faces(self.mesh.mesh, self.mesh.faces)
        if reverse_orientation:
            sgm.triangulate.reverse_face_orientations(self.mesh.mesh, self.mesh.faces)

        return self.mesh