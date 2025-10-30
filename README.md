# Seagull Mesh

Python bindings to CGAL's [Surface Mesh](https://doc.cgal.org/latest/Surface_mesh/index.html)
 and [Polygon Mesh Processing](https://doc.cgal.org/latest/Polygon_mesh_processing/index.html) modules, plus some other 
assorted extras.

## Installation

```shell
conda create --name seagull -c conda-forge cgal-cpp pybind11 eigen gmp mpfr 
conda activate seagull
```

Pyvista is used as a plotting utility in several bindings, but is completely optional. This can be installed
into the conda environment via pip:

```shell
pip install pyvista
```

Ceres-solver is used in the mesh smoothing module. This module functions without ceres-solver as well, so this dependency is optional. This can be installed with conda:

```shell
conda install -c conda-forge ceres-solver
```

On linux, you'll also need

```shell
conda install -c conda-forge gxx_linux-64 libgcc libxcrypt
```

Finally, 

```shell
git clone git@github.com:darikg/seagullmesh.git
pip install -vv ./seagullmesh
```

## Usage

### Overview

This package provides bindings to the `Surface_mesh<Point_3>` class using `Exact_predicates_inexact_constructions_kernel`. Because CGAL is a heavily templated library, only a few convenience bindings per overloaded method are available, more added as-needed. To simplify some of the python/C++ mappings, the bound classes and methods are further wrapped in a python layer.

### Mesh IO

3 constructors are available:
```python
from seagullmesh import Mesh3

mesh = Mesh3.from_polygon_soup(verts, faces)
mesh = Mesh3.from_file('mesh.ply')
mesh = Mesh3.from_pyvista(polydata)
```

with the corresponding output methods `mesh.to_polygon_soup`, `mesh.to_file`, `mesh.to_pyvista`.

#### Mesh construction from point clouds
Meshes can also be constructed from oriented point clouds using CGAL's [Poisson Surface Reconstruction package](
https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html)

```python
mesh = Mesh3.from_poisson_surface_reconstruction(points, normals, spacing)
```

Or from unoriented point clouds using CGAL's [3D Alpha Wrapping](https://doc.cgal.org/latest/Alpha_wrap_3/index.html).

```python
mesh = Mesh3.from_alpha_wrapping(points, alpha, offset)
```

### Mesh indices

CGAL indices `Surface_mesh::vertex_index`, `::edge_index`, `::face_index`, `::edge_index`, and `::halfedge_index` 
are exposed in the python properties `mesh.vertices`, `.faces`, `.edges`, and `.halfedges`, which are returned as 
opaque wrappers around `std::vectors` that allow numpy-like indexing and slicing. These vectors are used for indexing 
property maps and specifying regions for further processing. (See below.)

### Mesh property maps

Property maps are stored in the python properties `mesh.vertex_data`, 
`edge_data`, etc. They can be manually created by specifying a default value:

```python
mesh.vertex_data.create('my_property1', default=1)
mesh.edge_data.create('my_bool', default=False)
```

Because it may be necessary to distinguish between C++ signed and unsigned integer value types, an additional `signed` 
kwarg is supported for integer types:

```python
mesh.vertex_data.create('my_uint', default=0, signed=False)
```

Properties can also be created by assigned an array directly to a new property map name, 
where a default value of 0 is assumed:
```python
mesh.vertex_data['my_property2'] = np.arange(mesh.n_vertices)
```

Indexing can use the appropriate keys, which returns a numpy array of property values for the corresponding indices:
```python
my_property_vals = mesh.vertex_data['my_property2'][mesh.vertices]
```

Numpy-like indexing is also supported for convenience:
```python
first_10 = mesh.vertex_data['my_property2'][:10]
every_other = mesh.vertex_data['my_property2'][(np.arange(mesh.n_vertices) % 2) == 0]
```

Non-scalar valued properties are supported in the form of CGAL's `Point_2`, `Point_3`, `Vector_2`, `Vector_3`, exposed
as the python classes `seagullmesh.Point2`, `Point3`, `Vector2` and `Vector3`.

```python
from seagullmesh import Point2

mesh.vertex_data.create('uv_map', default=Point2(0, 0))
mesh.vertex_data['uv_map'] = np.random.uniform(-1, 1, (mesh.n_vertices, 2))
```

## Mesh processing

Currently implemented are:

From [PMP Meshing](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html):
  - `mesh.remesh(target_edge_length, ...)`
    - Isotropic remeshing
  - `mesh.remesh_adaptive(edge_len_min_max, tolerance, ...)`
    - Curvature-adaptive isotropic remeshing
  - `mesh.fair(vertices, fairing_continuity)`
  - `mesh.refine(faces, density)`
  - `mesh.smooth_angle_and_area(faces, n_iterations, do_angle_smoothing, do_area_smoothing)`
  - `mesh.smooth_shape(faces, time)`
  - `mesh.tangential_relaxation(verts)`

From [PMP Corrected Curvature Computation](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corrected__curvatures__grp.html)
  - `mesh.interpolated_corrected_curvatures(ball_radius, mc_map, gc_map, pc_map)`

From [PMP Corefinement and Boolean Operations](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html)
  - `mesh.corefine(other)`
  - `mesh.union(other)`
  - `mesh.difference(other)`
  - `mesh.intersection(other)`

From the [PMP Location Functions](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__locate__grp.html)
  - `mesh.aabb_tree()`
  - `mesh.locate_points(points, aabb_tree)`
  - `mesh.construct_points(faces, bary_coords)`

From [Triangulated Surface Mesh Shortest Paths
](https://doc.cgal.org/latest/Surface_mesh_shortest_path/group__PkgSurfaceMeshShortestPathRef.html)
  - `mesh.shortest_path(src_face, src_bary_coords, tgt_face, tgt_bary_coords)`

From [The Heat Method](https://doc.cgal.org/latest/Heat_method_3/classCGAL_1_1Heat__method__3_1_1Surface__mesh__geodesic__distances__3.html)
  - `mesh.estimate_geodesic_distances(source_vertex_or_vertices)`

From [Planar Parameterization of Triangulated Surface Meshes](https://doc.cgal.org/latest/Surface_mesh_parameterization/group__PkgSurfaceMeshParameterizationRef.html)
  - `mesh.lscm(uv_map)`
  - `mesh.arap(uv_map)`

From [Triangulated Surface Mesh Simplification](https://doc.cgal.org/latest/Surface_mesh_simplification/index.html)
  - `mesh.edge_collapse()`

From [Triangulated Surface Mesh Skeletonization](https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html)
  - `mesh.skeletonize()`

## Acknowledgements

The basis of the pybind11-cgal infrastructure is inspired very heavily by [Scikit-geometry](https://github.com/scikit-geometry/scikit-geometry).

## See also
  - [Scikit-geometry](https://github.com/scikit-geometry/scikit-geometry)
  - [pygalmesh](https://github.com/meshpro/pygalmesh)
  - [potpourri3d](https://github.com/nmwsharp/potpourri3d)
  - [pyvista](https://github.com/pyvista/pyvista)

## License

This software is licensed under the LGPL-3 license. See the [LICENSE](LICENSE) file for details.
