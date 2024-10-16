from pathlib import Path

from seagullmesh import Mesh3
from seagullmesh import _seagullmesh as sgm

file = Path(__file__).parent / 'test' / 'assets' / 'armadillo.off'


# m.interpolated_corrected_curvatures()
# print(m.vertex_data['principal_curvature'][m.vertices])

m = Mesh3.from_file(str(file))
edge_len = m.to_pyvista().extract_all_edges().compute_cell_sizes().cell_data['Length'].mean()

# m.remesh(
#     target_edge_length=edge_len * 2,
#     touched_map=None,
# )

m.remesh_adaptive(
    edge_len_min_max=(edge_len * 0.8, edge_len * 2),
    tolerance=edge_len / 10,
    n_iter=1,
    ball_radius=edge_len * 2,
    # touched_map='what',
)

# m.remesh()
# sizing = sgm.meshing.AdaptiveSizingField(tol=0.1, edge_len_min_max=(0.3, 0.5), faces=m.faces, mesh=m.mesh)