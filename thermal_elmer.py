#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""`pc cae thermal`, by Elmer: steady heat conduction through one part.

The `thermal:` keys, applied to the faces each port covers:

* `temperature:` holds them at a temperature (Celsius);
* `heat:` puts a power in through them (watts), spread evenly over their area;
* `convection:` cools them by the surrounding fluid with a film coefficient
  (W/(m^2*K)), against `ambient_temperature`.

Every other surface is insulated. A steady state needs somewhere for the heat
to go, so at least one surface has to be held at a temperature or cooled.

The answer is the temperature field, as the model coloured by it, and its range
as a measurement; a part hotter than `max_temperature` anywhere is a finding.
"""

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
    held, held_records, held_empty = support.condition_groups(request, "temperature", "temperature")
    heated, heated_records, heated_empty = support.condition_groups(request, "heat", "heat")
    cooled, cooled_records, cooled_empty = support.condition_groups(request, "convection", "convection")
    findings.extend(elmer.nothing_covered(record) for record in held_empty + heated_empty + cooled_empty)
    if not held_records and not cooled_records:
        raise support.ToolFailed(
            "'thermal:' gives the heat nowhere to go: hold a surface at a 'temperature:', or cool one by "
            "'convection:', on a port that covers a surface of the part. Without one there is no steady state"
        )

    with tempfile.TemporaryDirectory() as work:
        # A directory to keep the case in, for whoever is debugging one by hand.
        work = os.environ.get("PARTCAD_ELMER_CASE") or work
        mesh, boundaries = elmer.mesh_part(
            request, {**held, **heated, **cooled}, work, int(request.get("mesh_order", 2))
        )
        elmer.solve(work, _sif(request, boundaries, held, held_records, heated, heated_records, cooled, cooled_records))
        points, cells, fields = elmer.read_vtu(elmer.find_vtu(work))

    temperature = [value[0] for value in fields.get("temperature", [])]
    if not temperature:
        raise support.ToolFailed("ElmerSolver produced no temperatures")
    support.write_glb(path, points, elmer.surface(points, cells), temperature)

    hottest, coldest = max(temperature), min(temperature)
    limit = request.get("max_temperature")
    if limit is not None and hottest > float(limit):
        findings.append(
            {
                "severity": WARNING,
                "message": "The part reaches %.4g C, past the %.4g C it may." % (hottest, float(limit)),
                "where": "hottest point",
            }
        )
    put_in = sum(float(record["heat"]) for record in heated_records)
    return {
        "success": True,
        "findings": findings,
        "warnings": [
            "Temperature from %.4g C to %.4g C, %.4g W put in, over %d elements"
            % (coldest, hottest, put_in, mesh.element_count)
        ],
    }


def _sif(request, boundaries, held, held_records, heated, heated_records, cooled, cooled_records):
    """HeatSolve, the conductivity, and a condition per group.

    Temperatures stay in Celsius: conduction and convection are linear in the
    temperature, so the scale's zero does not matter -- only radiation would
    need kelvin, and nothing here radiates.
    """
    ambient = float(request.get("ambient_temperature", 20.0))
    sif = (
        elmer.header()
        + """
Equation 1
  Active Solvers(2) = 1 2
End
Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
%s
End
%s
Material 1
  Heat Conductivity = %.10g
  Density = %.10g
End
Body Force 1
End
"""
        % (
            elmer.SIF_LINEAR,
            elmer.result_output(2),
            float(request.get("conductivity", 50.0)),
            float(request.get("density", 7850)),
        )
    )
    number = 1
    for name, record in zip(held, held_records):
        sif += "Boundary Condition %d\n  %s\n  Temperature = %.10g\nEnd\n" % (
            number,
            elmer.targets(boundaries, [name]),
            float(record["temperature"]),
        )
        number += 1
    for name, record in zip(heated, heated_records):
        # A power over the area it goes in through: a flux, into the part.
        flux = float(record["heat"]) / (float(record["area"]) * elmer.M2_PER_MM2)
        sif += "Boundary Condition %d\n  %s\n  Heat Flux BC = True\n  Heat Flux = %.10g\nEnd\n" % (
            number,
            elmer.targets(boundaries, [name]),
            flux,
        )
        number += 1
    for name, record in zip(cooled, cooled_records):
        sif += (
            "Boundary Condition %d\n  %s\n  Heat Flux BC = True\n  Heat Transfer Coefficient = %.10g\n"
            "  External Temperature = %.10g\nEnd\n"
            % (number, elmer.targets(boundaries, [name]), float(record["convection"]), ambient)
        )
        number += 1
    return sif


# Importable as a module by its tests, and executed by PartCAD's meta-wrapper
# under the run name below -- so nothing here runs on import.
if __name__ == "__partcad_export__":
    pass
