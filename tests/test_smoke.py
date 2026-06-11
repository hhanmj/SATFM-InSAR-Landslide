import numpy as np

from satfm.manuscript_table2 import (
    build_scene_geometry,
    make_blind_scene,
    solve_one_seed,
)


def test_manuscript_solvers_return_finite_values_on_small_scene():
    scene = make_blind_scene(ny=24, nx=24)
    geom = build_scene_geometry(scene)
    estimates, diagnostics = solve_one_seed(scene, geom, seed=20260424)
    n_pix = 24 * 24
    for method in ["WLS", "APFM", "SPFM", "Full SATFM"]:
        estimate = estimates[method]
        assert estimate.shape == (n_pix, 3)
        assert np.isfinite(estimate).all()
    assert diagnostics["edge_count"] > 0
