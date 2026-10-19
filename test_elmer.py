#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""What can be checked about this package, with and without Elmer in the room.

Without it: the VTU reader (ASCII parsing of a file Elmer does not spell
consistently), the eigenvalue reader, the solver input files each analysis
writes, and that a machine with no Elmer is told what to do.

With it -- in this package's image, or on a workstation with Elmer, gmsh and
build123d -- the example cantilever and pipe are solved and compared with their
closed forms: the analyses themselves, end to end.

Run with `pytest` from this directory, beside an installed `partcad` (for the
sandbox helpers the scripts import; see `conftest.py`).
"""

import math
import os

import pytest

support = pytest.importorskip("cae_support", reason="needs a partcad with wrappers/cae_support.py")

import cfd_elmer  # noqa: E402
import dynamic_elmer  # noqa: E402
import elmer_common as elmer  # noqa: E402
import fea_elmer  # noqa: E402
import thermal_elmer  # noqa: E402

VTU = """<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="5" NumberOfCells="2">
      <PointData>
        <DataArray type="Float64" Name="temperature" NumberOfComponents="1" format="ascii">
 1 2 3 4 5
        </DataArray>
        <DataArray type="Float64" Name="displacement EigenMode1" NumberOfComponents="3" format="ascii">
 0 0 0  0 0 1  0 0 2  0 0 3  0 0 4
        </DataArray>
      </PointData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
 0 0 0  1 0 0  0 1 0  0 0 1  1 1 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
 0 1 2 3  1 2 3 4
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
 4 8
        </DataArray>
        <DataArray type="Int32" Name="types" FORMAT="ascii">
 10 10
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"""


def test_the_vtu_reader_reads_what_elmer_writes(tmp_path):
    path = tmp_path / "result_t0001.vtu"
    path.write_text(VTU)
    assert elmer.find_vtu(str(tmp_path)) == str(path)
    points, cells, fields = elmer.read_vtu(str(path))
    assert points[4] == (1.0, 1.0, 1.0)
    assert cells == {10: [[0, 1, 2, 3], [1, 2, 3, 4]]}
    assert [value[0] for value in fields["temperature"]] == [1, 2, 3, 4, 5]
    # The shared face of the two tetrahedra is not part of the surface.
    assert len(elmer.surface(points, cells)) == 6


def test_a_binary_vtu_is_refused_with_a_reason(tmp_path):
    path = tmp_path / "result.vtu"
    path.write_text(VTU.replace('format="ascii"', 'format="appended"'))
    with pytest.raises(support.ToolFailed, match="binary"):
        elmer.read_vtu(str(path))


def test_the_first_mode_is_scaled_to_one():
    fields = {"displacement EigenMode1": [(0, 0, 0), (0, 0, 2), (0, 0, 4)]}
    assert dynamic_elmer._first_mode(fields) == [0.0, 0.5, 1.0]


def test_eigenvalues_are_read_from_the_log():
    log = """
 EigenSolve: Computed 3 Eigen Values
 EigenSolve: --------------------------------
 EigenSolve: 1:   2.6109E+07   0.0000E+00
 EigenSolve: 2:   2.6110E+07   0.0000E+00
 EigenSolve: 3:   9.3929E+08   0.0000E+00
"""
    assert elmer.eigenvalues(log) == [2.6109e7, 2.6110e7, 9.3929e8]


def test_a_machine_with_no_elmer_is_told_what_to_install(monkeypatch):
    monkeypatch.setattr(elmer, "_found", None)
    monkeypatch.delenv(elmer.BIN_ENV, raising=False)
    monkeypatch.delenv("ELMER_HOME", raising=False)
    monkeypatch.setattr(elmer.shutil, "which", lambda _name: None)
    monkeypatch.setattr(elmer, "DIRECTORIES", ())
    with pytest.raises(support.ToolMissing) as raised:
        elmer.find_elmer()
    message = str(raised.value)
    assert "dockerImage:" in message
    assert "elmerfem-csc" in message
    assert elmer.BIN_ENV in message


def test_a_load_is_a_traction_along_minus_z():
    sif = fea_elmer._sif(
        {"self_weight": False},
        {"fix0": 1, "load0": 2},
        {"fix0": [4]},
        {"load0": [5]},
        [{"load": 100.0, "area": 100.0}],
    )
    # 100 N over 100 mm^2 is 1 MPa.
    assert "Force 3 = -1000000" in sif
    assert "Displacement 1 = 0" in sif
    assert "Stress Bodyforce 3 = 0\n" in sif


def test_a_load_with_a_direction_is_a_traction_along_it():
    sif = fea_elmer._sif(
        {"self_weight": False},
        {"fix0": 1, "load0": 2},
        {"fix0": [4]},
        {"load0": [5]},
        [{"load": 5.0, "force": [3.0, 0.0, -4.0], "area": 1.0}],
    )
    # 3 N and -4 N over 1 mm^2; nothing along Y, so no Y at all.
    assert "Force 1 = 3000000" in sif and "Force 3 = -4000000" in sif
    assert "Force 2" not in sif


def test_the_part_is_weighed_under_the_gravity_it_states():
    sif = fea_elmer._sif({"density": 1000, "gravity": [0.0, -8.0, -6.0]}, {"fix0": 1}, {"fix0": [4]}, {}, [])
    assert "Stress Bodyforce 1 = 0" in sif
    assert "Stress Bodyforce 2 = -8000" in sif and "Stress Bodyforce 3 = -6000" in sif
    # Its own gravity replaces the default, and nought weighs nothing.
    sif = fea_elmer._sif({"density": 1000, "gravity": [0, 0, 0]}, {"fix0": 1}, {"fix0": [4]}, {}, [])
    assert "Stress Bodyforce 3 = 0" in sif
    sif = fea_elmer._sif({"density": 1000}, {"fix0": 1}, {"fix0": [4]}, {}, [])
    assert "Stress Bodyforce 3 = -9810" in sif


def test_the_fea_file_type_says_it_reads_a_direction_and_a_gravity():
    import yaml

    with open(os.path.join(os.path.dirname(__file__), "partcad.yaml"), encoding="utf-8") as f:
        declared = yaml.safe_load(f)["cae"]["fea"]
    assert sorted(declared["reads"]) == ["force", "gravity"]


def test_heat_is_a_flux_and_convection_meets_the_ambient():
    sif = thermal_elmer._sif(
        {"ambient_temperature": 25.0},
        {"temperature0": 1, "heat0": 2, "convection0": 3},
        {"temperature0": [1]},
        [{"temperature": 20.0}],
        {"heat0": [2]},
        [{"heat": 1.0, "area": 100.0}],
        {"convection0": [3]},
        [{"convection": 10.0}],
    )
    assert "Temperature = 20" in sif
    assert "Heat Flux = 10000" in sif
    assert "Heat Transfer Coefficient = 10" in sif and "External Temperature = 25" in sif


def test_an_inlet_pushes_in():
    sif = cfd_elmer._sif({}, {"inlet0": 1, "outlet0": 2, "walls": 3}, ["inlet0"], [63.66], ["outlet0"], 900.0, 0.1)
    # A traction along the outward normal: pushing in is negative.
    assert "External Pressure = -63.66" in sif
    assert "Noslip Wall BC = True" in sif


def test_nothing_held_is_said_rather_than_solved():
    for module in (fea_elmer, dynamic_elmer):
        result = module.process("/nonexistent.glb", {"boundary": []})
        assert result["success"] is False and "holds nothing" in result["exception"]
    result = thermal_elmer.process("/nonexistent.glb", {"boundary": [{"heat": 1.0, "faces": [{"index": 0}]}]})
    assert result["success"] is False and "nowhere to go" in result["exception"]


# --------------------------------------------------------------------------- #
# The analyses themselves, where Elmer is                                     #
# --------------------------------------------------------------------------- #


def _elmer_here():
    try:
        elmer.find_elmer()
        return True
    except support.ToolMissing:
        return False


needs_elmer = pytest.mark.skipif(not _elmer_here(), reason="needs Elmer")


def _cantilever(analysis):
    build123d = pytest.importorskip("build123d")
    pytest.importorskip("gmsh")
    root = {"port": "root", "location": [[-50, 0, 0], [0, 0, 1], 0]}
    tip = {"port": "tip", "location": [[50, 0, 0], [0, 0, 1], 0]}
    if analysis == "fea":
        boundary = [dict(root, fix=True, load=0.0), dict(tip, fix=False, load=100.0)]
    elif analysis == "dynamic":
        boundary = [dict(root, fix=True)]
    else:
        boundary = [
            dict(root, temperature=20.0, heat=None, convection=None),
            dict(tip, temperature=None, heat=1.0, convection=None),
        ]
    request = {
        "analysis": analysis,
        "wrapped": build123d.Box(100, 10, 10).wrapped,
        "boundary": boundary,
        "youngs_modulus": 2.0e11,
        "poissons_ratio": 0.3,
        "density": 7850,
        "self_weight": False,
        "mesh_size": 0.03,
        "mesh_order": 2,
        "conductivity": 50.0,
        "modes": 2,
    }
    assert support.prepare(request) == []
    return request


def _measured(result, label):
    assert result["success"], result.get("exception")
    return result["warnings"][0].split(label)[1].split()[0].rstrip(",")


@needs_elmer
def test_the_cantilever_deflects_as_a_beam(tmp_path):
    """100 N on a 100x10x10 mm steel bar: PL^3/3EI = 0.2000 mm."""
    result = fea_elmer.process(str(tmp_path / "c.glb"), _cantilever("fea"))
    assert float(_measured(result, "peak displacement")) == pytest.approx(0.2016, rel=0.01)


@needs_elmer
def test_the_cantilever_rings_as_a_beam(tmp_path):
    """The first bending frequency: (1.875^2 / 2 pi L^2) sqrt(EI / rho A) = 815.4 Hz."""
    result = dynamic_elmer.process(str(tmp_path / "c.glb"), _cantilever("dynamic"))
    first = float(_measured(result, "Natural frequencies"))
    expected = 1.875104**2 / (2 * math.pi * 0.1**2) * math.sqrt(2.0e11 * 0.01**4 / 12 / (7850 * 1e-4))
    assert first == pytest.approx(expected, rel=0.01)


@needs_elmer
def test_the_cantilever_warms_as_a_rod(tmp_path):
    """1 W through a 10x10 mm steel bar 100 mm long: P L / k A = 20 K over the root."""
    result = thermal_elmer.process(str(tmp_path / "c.glb"), _cantilever("thermal"))
    assert result["success"], result.get("exception")
    hottest = float(result["warnings"][0].split(" C to ")[1].split(" C")[0])
    assert hottest == pytest.approx(40.0, rel=1e-3)
    assert os.path.getsize(tmp_path / "c.glb") > 0
