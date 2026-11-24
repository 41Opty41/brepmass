"""Tests for the Python ctypes wrapper around libbrepmass."""

from __future__ import annotations

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "python"))
import brepmass


def test_python_tessellate_box_and_mass_props() -> None:
    verts, tris = brepmass.tessellate_box(2.0, 3.0, 4.0)
    props = brepmass.mesh_mass_props(verts, tris, density=2.0)
    assert props.dim == 3
    assert props.is_closed is True
    assert props.euler_char == 2
    assert props.volume == pytest.approx(24.0, rel=1e-12)
    assert props.surface_area == pytest.approx(52.0, rel=1e-12)
    assert props.mass == pytest.approx(48.0, rel=1e-12)
    assert props.centroid == pytest.approx((1.0, 1.5, 2.0), abs=1e-12)
    # Principal moments of 2x3x4 box of mass 48: m*(b^2+c^2)/12
    assert props.principal_moments == pytest.approx(
        (48.0 * (4.0 + 9.0) / 12.0, 48.0 * (4.0 + 16.0) / 12.0, 48.0 * (9.0 + 16.0) / 12.0),
        rel=1e-12,
    )
