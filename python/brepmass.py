"""Thin Python ctypes wrapper for the C99 brepmass shared library."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import os
from pathlib import Path
import sys
from typing import Sequence


class _CMassProps(ctypes.Structure):
    _fields_ = [
        ("dim", ctypes.c_int),
        ("is_closed", ctypes.c_int),
        ("euler_char", ctypes.c_int),
        ("n_vertices", ctypes.c_int),
        ("n_triangles", ctypes.c_int),
        ("volume", ctypes.c_double),
        ("surface_area", ctypes.c_double),
        ("mass", ctypes.c_double),
        ("centroid", ctypes.c_double * 3),
        ("inertia_origin", (ctypes.c_double * 3) * 3),
        ("inertia_centroid", (ctypes.c_double * 3) * 3),
        ("principal_moments", ctypes.c_double * 3),
        ("principal_axes", (ctypes.c_double * 3) * 3),
        ("bbox_min", ctypes.c_double * 3),
        ("bbox_max", ctypes.c_double * 3),
    ]


class _CMesh(ctypes.Structure):
    _fields_ = [
        ("vertices", ctypes.POINTER(ctypes.c_double)),
        ("triangles", ctypes.POINTER(ctypes.c_int)),
        ("face_ids", ctypes.POINTER(ctypes.c_int)),
        ("n_vertices", ctypes.c_int),
        ("cap_vertices", ctypes.c_int),
        ("n_triangles", ctypes.c_int),
        ("cap_triangles", ctypes.c_int),
    ]


class _CTessOpts(ctypes.Structure):
    _fields_ = [
        ("u_segments", ctypes.c_int),
        ("v_segments", ctypes.c_int),
        ("chord_tol", ctypes.c_double),
    ]


@dataclass
class MassProps:
    """Computed mass, centroid, inertia tensors, and topology report."""

    dim: int
    is_closed: bool
    euler_char: int
    n_vertices: int
    n_triangles: int
    volume: float
    surface_area: float
    mass: float
    centroid: tuple[float, float, float]
    inertia_origin: tuple[tuple[float, float, float], ...]
    inertia_centroid: tuple[tuple[float, float, float], ...]
    principal_moments: tuple[float, float, float]
    principal_axes: tuple[tuple[float, float, float], ...]
    bbox_min: tuple[float, float, float]
    bbox_max: tuple[float, float, float]


def _find_library() -> Path:
    env = os.environ.get("BREPMASS_LIB")
    if env and Path(env).exists():
        return Path(env)
    root = Path(__file__).resolve().parent.parent
    names = ["brepmass.dll", "libbrepmass.so", "libbrepmass.dylib"]
    for cand_dir in (root / "build", root / "build" / "Release", root):
        for nm in names:
            p = cand_dir / nm
            if p.exists():
                return p
    raise FileNotFoundError(
        "brepmass shared library not found; build with CMake/clang or set BREPMASS_LIB."
    )


def _load_lib() -> ctypes.CDLL:
    lib = ctypes.CDLL(str(_find_library()))
    lib.bm_mesh_mass_props.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.POINTER(_CMassProps),
    ]
    lib.bm_mesh_mass_props.restype = ctypes.c_int

    lib.bm_mesh_init.argtypes = [ctypes.POINTER(_CMesh)]
    lib.bm_mesh_init.restype = None
    lib.bm_mesh_free.argtypes = [ctypes.POINTER(_CMesh)]
    lib.bm_mesh_free.restype = None

    lib.bm_tessellate_box.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(_CMesh)]
    lib.bm_tessellate_box.restype = ctypes.c_int

    lib.bm_tessellate_sphere.argtypes = [ctypes.c_double, ctypes.POINTER(_CTessOpts), ctypes.POINTER(_CMesh)]
    lib.bm_tessellate_sphere.restype = ctypes.c_int
    return lib


def _to_mass_props(raw: _CMassProps) -> MassProps:
    return MassProps(
        dim=raw.dim,
        is_closed=bool(raw.is_closed),
        euler_char=raw.euler_char,
        n_vertices=raw.n_vertices,
        n_triangles=raw.n_triangles,
        volume=raw.volume,
        surface_area=raw.surface_area,
        mass=raw.mass,
        centroid=(raw.centroid[0], raw.centroid[1], raw.centroid[2]),
        inertia_origin=tuple((row[0], row[1], row[2]) for row in raw.inertia_origin),
        inertia_centroid=tuple((row[0], row[1], row[2]) for row in raw.inertia_centroid),
        principal_moments=(raw.principal_moments[0], raw.principal_moments[1], raw.principal_moments[2]),
        principal_axes=tuple((row[0], row[1], row[2]) for row in raw.principal_axes),
        bbox_min=(raw.bbox_min[0], raw.bbox_min[1], raw.bbox_min[2]),
        bbox_max=(raw.bbox_max[0], raw.bbox_max[1], raw.bbox_max[2]),
    )


def mesh_mass_props(
    vertices: Sequence[Sequence[float]],
    triangles: Sequence[Sequence[int]],
    *,
    density: float = 1.0,
) -> MassProps:
    """Compute exact volume, surface area, centroid, and inertia tensors via C99 brepmass."""
    lib = _load_lib()
    nv = len(vertices)
    nt = len(triangles)
    v_flat = (ctypes.c_double * (3 * nv))(*(float(c) for p in vertices for c in p[:3]))
    t_flat = (ctypes.c_int * (3 * nt))(*(int(i) for tri in triangles for i in tri[:3]))
    out = _CMassProps()
    rc = lib.bm_mesh_mass_props(v_flat, nv, t_flat, nt, float(density), ctypes.byref(out))
    if rc != 0:
        raise RuntimeError(f"bm_mesh_mass_props returned error code {rc}")
    return _to_mass_props(out)


def tessellate_box(dx: float, dy: float, dz: float) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    """Tessellate a box via C99 brepmass and return (vertices, triangles)."""
    lib = _load_lib()
    cm = _CMesh()
    lib.bm_mesh_init(ctypes.byref(cm))
    try:
        rc = lib.bm_tessellate_box(float(dx), float(dy), float(dz), ctypes.byref(cm))
        if rc != 0:
            raise RuntimeError(f"bm_tessellate_box failed ({rc})")
        verts = [
            (cm.vertices[3 * i], cm.vertices[3 * i + 1], cm.vertices[3 * i + 2])
            for i in range(cm.n_vertices)
        ]
        tris = [
            (cm.triangles[3 * k], cm.triangles[3 * k + 1], cm.triangles[3 * k + 2])
            for k in range(cm.n_triangles)
        ]
        return verts, tris
    finally:
        lib.bm_mesh_free(ctypes.byref(cm))
