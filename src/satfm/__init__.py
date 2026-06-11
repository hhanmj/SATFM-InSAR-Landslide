"""SATFM manuscript Table 2 synthetic benchmark package."""

from .manuscript_table2 import (
    build_edges,
    build_scene_geometry,
    make_blind_scene,
    run_manuscript_table2,
    solve_adaptive_global,
    solve_strict_apfm,
    solve_strict_spfm,
    solve_weighted_local,
)

__all__ = [
    "build_edges",
    "build_scene_geometry",
    "make_blind_scene",
    "run_manuscript_table2",
    "solve_adaptive_global",
    "solve_strict_apfm",
    "solve_strict_spfm",
    "solve_weighted_local",
]
