#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae fea`, by Elmer: a linear static stress analysis of one part.

`fix:` holds the faces its ports cover still in every direction; `load:` spreads
its force evenly over the faces its ports cover, along -Z, which is where a
weight goes -- PartCAD's `load:` is a magnitude, and a part cannot state a
direction yet. Gravity acts on the part's own mass too, unless `self_weight` is
off.

Quadratic tetrahedra by default: linear ones are far too stiff in bending, which
is the loading almost every bracket sees.
"""

import os
import sys
import tempfile

# `runpy.run_path()` leaves the script's directory off `sys.path`, so the
# sibling module below is not importable without saying where it is.
sys.path.append(os.path.dirname(__file__))
import elmer_common as elmer  # noqa: E402
from elmer_common import ERROR, WARNING, support  # noqa: E402

MM_PER_M = 1000.0


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    fixed, fixed_records, fixed_empty = support.condition_groups(request, "fix", "fix")
    loaded, loaded_records, loaded_empty = support.condition_groups(request, "load", "load")
    findings.extend(elmer.nothing_covered(record) for record in fixed_empty + loaded_empty)
    if not fixed_records:
        # Nothing to react against: every node would move together and the
        # stiffness matrix is singular.
        raise support.ToolFailed(
            "'fea:' holds nothing, so the part is free to float: name at least one interface under 'fix:', "
            "on a port that covers a surface of the part"
        )

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_ELMER_CASE") or work
        mesh, boundaries = elmer.mesh_part(request, {**fixed, **loaded}, work, int(request.get("mesh_order", 2)))
        elmer.solve(work, _sif(request, boundaries, fixed, loaded, loaded_records))
        points, cells, fields = elmer.read_vtu(elmer.find_vtu(work))

    displacement = elmer.magnitude(fields.get("displacement", []))
    stress = [value[0] for value in fields.get("vonmises", [])]
    if not displacement or not stress:
        raise support.ToolFailed("ElmerSolver produced no displacements or no stresses")
    support.write_glb(path, points, elmer.surface(points, cells), stress, low=0.0)

    peak_stress = max(stress)
    peak_displacement = max(displacement) * MM_PER_M
    findings.extend(_verdict(peak_stress, peak_displacement, request))
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Peak von Mises stress %.4g MPa, peak displacement %.4g mm, over %d elements"
            % (peak_stress / 1e6, peak_displacement, mesh.element_count)
        ],
    }


def _sif(request, boundaries, fixed, loaded, loaded_records):
    """The solver input file: StressSolve, the material, gravity, and a condition per group."""
    density = float(request.get("density", 7850))
    gravity = 9.81 if request.get("self_weight", True) else 0.0
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
  Calculate Stresses = True
%s
End
%s
Material 1
  Youngs Modulus = %.10g
  Poisson Ratio = %.10g
  Density = %.10g
End
Body Force 1
  Stress Bodyforce 3 = %.10g
End
"""
        % (
            elmer.SIF_LINEAR,
            elmer.result_output(2),
            float(request.get("youngs_modulus", 2.1e11)),
            float(request.get("poissons_ratio", 0.30)),
            density,
            -density * gravity,
        )
    )
    number = 1
    for name in fixed:
        sif += (
            "Boundary Condition %d\n  %s\n  Displacement 1 = 0\n  Displacement 2 = 0\n  Displacement 3 = 0\nEnd\n"
            % (
                number,
                elmer.targets(boundaries, [name]),
            )
        )
        number += 1
    for name, record in zip(loaded, loaded_records):
        # A force on the port, spread over the area it covers: a traction, so a
        # finer mesh does not load the part harder.
        traction = float(record["load"]) / (float(record["area"]) * elmer.M2_PER_MM2)
        sif += "Boundary Condition %d\n  %s\n  Force 3 = %.10g\nEnd\n" % (
            number,
            elmer.targets(boundaries, [name]),
            -traction,
        )
        number += 1
    return sif


def _verdict(peak_stress, peak_displacement, request):
    """What the numbers mean, as findings: the thresholds the package configured."""
    findings = []
    yield_strength = float(request.get("yield_strength", 2.5e8))
    safety_factor = float(request.get("safety_factor", 1.5))
    allowable = yield_strength / max(safety_factor, 1e-9)
    if peak_stress > yield_strength:
        findings.append(
            {
                "severity": ERROR,
                "message": "The part yields under this load: peak stress %.3g MPa against a yield strength of %.3g MPa."
                % (peak_stress / 1e6, yield_strength / 1e6),
                "where": "peak von Mises stress",
            }
        )
    elif peak_stress > allowable:
        findings.append(
            {
                "severity": WARNING,
                "message": "The part does not reach the safety factor of %.2f: peak stress %.3g MPa against an "
                "allowable %.3g MPa." % (safety_factor, peak_stress / 1e6, allowable / 1e6),
                "where": "peak von Mises stress",
            }
        )
    limit = request.get("max_displacement")
    if limit is not None and peak_displacement > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part moves %.3g mm under this load, past the %.3g mm it is allowed."
                % (peak_displacement, float(limit)),
                "where": "peak displacement",
            }
        )
    return findings


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
