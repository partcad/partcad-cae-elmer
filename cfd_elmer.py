#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae cfd`, by Elmer: steady, laminar, incompressible flow through one part.

The part is read as the **fluid volume**. Its `cfd:` keys, applied to the faces
each port covers:

* `load:` names an **inlet**, pushed on with the pressure that is its force over
  the area its port covers;
* `outlet:` names where the flow **leaves**, at the reference pressure;
* every other surface is a **wall**, no-slip.

At inlets and outlets the flow is held normal to the surface (no tangential
velocity), which is what a developed flow through a cut does and what keeps a
pressure boundary from being met by a swirl it cannot answer.

FlowSolve, stabilized linear elements, Newton after the first iterations. The
flow rate is integrated over the outlets from the velocities on them; the model
is the boundary coloured by pressure.
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

TURBULENT_REYNOLDS = 2300.0


def process(path, request):
    try:
        return _analyse(path, request)
    except (support.ToolMissing, support.ToolFailed) as e:
        return support.failed(e)


def _analyse(path, request):
    findings = []
    density = float(request.get("density", 1.204))
    viscosity = float(request.get("viscosity", 1.82e-5))
    inlets, inlet_records, inlet_empty = support.condition_groups(request, "load", "inlet")
    outlets, outlet_records, outlet_empty = support.condition_groups(request, "outlet", "outlet")
    findings.extend(elmer.nothing_covered(record) for record in inlet_empty + outlet_empty)
    if not inlet_records:
        raise support.ToolFailed("'cfd:' drives nothing: name an inlet under 'load:', on a port that covers a surface")
    if not outlet_records:
        raise support.ToolFailed(
            "'cfd:' names no outlet, so the flow has nowhere to go and the problem has no downstream pressure "
            "reference: name where it leaves under 'outlet:', on a port that covers a surface"
        )

    pressures = [float(r["load"]) / (float(r["area"]) * elmer.M2_PER_MM2) for r in inlet_records]
    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_ELMER_CASE") or work
        mesh, boundaries = elmer.mesh_part(request, {**inlets, **outlets}, work, 1, rest="walls")
        log = elmer.solve(work, _sif(request, boundaries, list(inlets), pressures, list(outlets), density, viscosity))
        points, cells, fields = elmer.read_vtu(elmer.find_vtu(work))

    velocity = fields.get("velocity")
    pressure = [value[0] for value in fields.get("pressure", [])]
    if not velocity or not pressure:
        raise support.ToolFailed("ElmerSolver produced no flow field")
    support.write_glb(path, points, elmer.surface(points, cells), pressure)

    # The triangles are gmsh's and the velocities Elmer's, joined by position.
    at = support.match_points(mesh.coordinates, points, mesh.extent * 1e-6)
    flow = sum(
        support.flow_through(
            mesh,
            name,
            lambda k: velocity[at[k]] if at[k] is not None else None,
            support.outward_normal(record),
        )
        for name, record in zip(outlets, outlet_records)
    )
    outlet_area = sum(float(record["area"]) for record in outlet_records) * elmer.M2_PER_MM2
    mean_speed = flow / outlet_area if outlet_area else 0.0
    peak_speed = max(elmer.magnitude(velocity))
    reynolds = density * abs(mean_speed) * math.sqrt(4.0 * outlet_area / math.pi) / viscosity

    if "did not converge" in log.lower() or "max iterations reached" in log.lower():
        findings.append(
            {"severity": WARNING, "message": "The flow did not settle: these numbers are not a steady state."}
        )
    if reynolds > TURBULENT_REYNOLDS:
        findings.append(
            {
                "severity": WARNING,
                "message": "The flow is likely turbulent (Reynolds number about %.0f, past %.0f), and this analysis "
                "is laminar." % (reynolds, TURBULENT_REYNOLDS),
            }
        )
    limit = request.get("max_velocity")
    if limit is not None and peak_speed > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The flow reaches %.3g m/s, past the %.3g m/s it is allowed." % (peak_speed, float(limit)),
                "where": "peak speed",
            }
        )
    minimum = request.get("min_flow_rate")
    if minimum is not None and flow < float(minimum):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part passes %.4g m^3/s, short of the %.4g m^3/s it has to." % (flow, float(minimum)),
                "where": "flow rate",
            }
        )
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Flow rate %.4g m^3/s (%.4g L/min), mean outlet speed %.4g m/s, peak speed %.4g m/s, Reynolds number "
            "%.3g, driven by %s, over %d elements"
            % (
                flow,
                flow * 60000.0,
                mean_speed,
                peak_speed,
                reynolds,
                ", ".join("%.4g Pa" % value for value in pressures),
                mesh.element_count,
            )
        ],
    }


def _sif(request, boundaries, inlets, pressures, outlets, density, viscosity):
    """FlowSolve with the inlets pushed on, the outlets free, and every wall held."""
    sif = (
        elmer.header()
        + """
Equation 1
  Active Solvers(2) = 1 2
End
Solver 1
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Variable = Flow Solution[Velocity:3 Pressure:1]
  Stabilization Method = Stabilized
  Nonlinear System Max Iterations = %d
  Nonlinear System Convergence Tolerance = 1.0e-8
  Nonlinear System Newton After Iterations = 2
  Nonlinear System Newton After Tolerance = 1.0e-4
%s
End
%s
Material 1
  Density = %.10g
  Viscosity = %.10g
End
Body Force 1
End
"""
        % (
            int(request.get("iterations", 50)),
            elmer.SIF_ITERATIVE,
            elmer.result_output(2),
            density,
            viscosity,
        )
    )
    number = 1
    for name, pressure in zip(inlets, pressures):
        # 'External Pressure' is a traction along the outward normal, so a
        # pressure pushing *in* is a negative one.
        sif += (
            "Boundary Condition %d\n  %s\n  Normal-Tangential Velocity = True\n  Velocity 2 = 0\n  Velocity 3 = 0\n"
            "  External Pressure = %.10g\nEnd\n" % (number, elmer.targets(boundaries, [name]), -pressure)
        )
        number += 1
    for name in outlets:
        sif += (
            "Boundary Condition %d\n  %s\n  Normal-Tangential Velocity = True\n  Velocity 2 = 0\n  Velocity 3 = 0\n"
            "End\n" % (number, elmer.targets(boundaries, [name]))
        )
        number += 1
    sif += "Boundary Condition %d\n  %s\n  Noslip Wall BC = True\nEnd\n" % (
        number,
        elmer.targets(boundaries, ["walls"]),
    )
    return sif


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
