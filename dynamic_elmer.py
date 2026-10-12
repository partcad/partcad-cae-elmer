#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae dynamic`, by Elmer: the natural frequencies of one part, held as it is.

A modal analysis: `fix:` holds the faces its ports cover still, and ElmerSolver
finds the lowest `modes` frequencies at which what is left rings on its own,
with the shape of each. Nothing is loaded -- a natural frequency is a property
of the part, its material and how it is held, not of what it carries.

The answer is the frequencies, as a measurement, and the shape of the first
mode as the model: a part coloured by how far each point moves when it rings at
its lowest frequency, the thing to stiffen if that frequency is too low. A
first frequency under `min_frequency` is a finding: a part that rings within
reach of whatever drives it -- a motor, a fan, a road -- resonates.
"""

import math
import os
import sys
import tempfile

# `runpy.run_path()` leaves the script's directory off `sys.path`, so the
# sibling module below is not importable without saying where it is.
sys.path.append(os.path.dirname(__file__))
import elmer_common as elmer  # noqa: E402
from elmer_common import WARNING, support  # noqa: E402


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    fixed, fixed_records, fixed_empty = support.condition_groups(request, "fix", "fix")
    findings.extend(elmer.nothing_covered(record) for record in fixed_empty)
    if not fixed_records:
        # A part held by nothing rings at zero: six rigid-body modes, and the
        # answer would be those. Said here rather than reported as six zeros.
        raise support.ToolFailed(
            "'dynamic:' holds nothing, so the part is free and its lowest modes are rigid-body motion at 0 Hz: "
            "name what holds it under 'fix:', on a port that covers a surface of the part"
        )
    modes = max(1, int(request.get("modes", 6)))

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_ELMER_CASE") or work
        mesh, boundaries = elmer.mesh_part(request, fixed, work, int(request.get("mesh_order", 2)))
        log = elmer.solve(work, _sif(request, boundaries, fixed, modes))
        points, cells, fields = elmer.read_vtu(elmer.find_vtu(work))

    eigenvalues = elmer.eigenvalues(log)[:modes]
    if not eigenvalues:
        raise support.ToolFailed("ElmerSolver computed no eigenvalues: the modal analysis did not converge")
    frequencies = [math.sqrt(max(value, 0.0)) / (2.0 * math.pi) for value in eigenvalues]

    shape = _first_mode(fields)
    support.write_glb(path, points, elmer.surface(points, cells), shape or [0.0] * len(points), low=0.0)

    limit = request.get("min_frequency")
    if limit is not None and frequencies[0] < float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part's first natural frequency is %.4g Hz, under the %.4g Hz it has to stay above: "
                "it will resonate with anything driving it near that." % (frequencies[0], float(limit)),
                "where": "first natural frequency",
            }
        )
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Natural frequencies %s Hz, over %d elements"
            % (", ".join("%.5g" % value for value in frequencies), mesh.element_count)
        ],
    }


def _first_mode(fields):
    """How far each point moves in the first mode, scaled so the most is 1."""
    for name in sorted(fields):
        if "eigenmode1" in name.lower().replace(" ", "") or name.lower().endswith("eigenmode1"):
            moved = elmer.magnitude(fields[name])
            peak = max(moved, default=0.0)
            return [value / peak for value in moved] if peak > 0 else moved
    return None


def _sif(request, boundaries, fixed, modes):
    """StressSolve as an eigen analysis: the lowest 'modes' modes, shift-and-invert about zero."""
    sif = (
        elmer.header()
        + """
Equation 1
  Active Solvers(2) = 1 2
End
Solver 1
  Equation = Linear elasticity
  Procedure = "StressSolve" "StressSolver"
  Variable = -dofs 3 Displacement
  Eigen Analysis = True
  Eigen System Values = %d
  Eigen System Select = Smallest Magnitude
  Eigen System Convergence Tolerance = 1.0e-9
%s
End
%s
Material 1
  Youngs Modulus = %.10g
  Poisson Ratio = %.10g
  Density = %.10g
End
Body Force 1
End
"""
        % (
            modes,
            elmer.SIF_LINEAR,
            elmer.result_output(2),
            float(request.get("youngs_modulus", 2.1e11)),
            float(request.get("poissons_ratio", 0.30)),
            float(request.get("density", 7850)),
        )
    )
    for number, name in enumerate(fixed, start=1):
        sif += (
            "Boundary Condition %d\n  %s\n  Displacement 1 = 0\n  Displacement 2 = 0\n  Displacement 3 = 0\nEnd\n"
            % (
                number,
                elmer.targets(boundaries, [name]),
            )
        )
    return sif


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
