#!/usr/bin/env python3
#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""Whether this image can actually run the Elmer analyses.

Run at build time, so that an image which cannot is never published. The base
image's own `verify.py` checks the contract PartCAD talks to a container
through; this checks what this image adds, the way the analyses use it: gmsh
meshes a bar from *this* interpreter, ElmerGrid converts the mesh, and
ElmerSolver solves heat conduction along it -- which also proves ElmerSolver
finds its own solver modules through ELMER_HOME.
"""

import os
import re
import subprocess
import sys
import tempfile


def main() -> int:
    failures = []
    with tempfile.TemporaryDirectory() as work:
        msh = os.path.join(work, "bar.msh")
        try:
            import gmsh

            gmsh.initialize()
            try:
                gmsh.option.setNumber("General.Terminal", 0)
                gmsh.model.add("verify")
                gmsh.model.occ.addBox(0, 0, 0, 1, 0.1, 0.1)
                gmsh.model.occ.synchronize()
                for dim, tag in gmsh.model.getEntities(2):
                    x = gmsh.model.occ.getCenterOfMass(dim, tag)[0]
                    if abs(x) < 1e-9:
                        gmsh.model.addPhysicalGroup(2, [tag], name="cold")
                    elif abs(x - 1) < 1e-9:
                        gmsh.model.addPhysicalGroup(2, [tag], name="hot")
                gmsh.model.addPhysicalGroup(3, [1], name="body")
                gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
                gmsh.model.mesh.generate(3)
                gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
                gmsh.write(msh)
                print("gmsh %s meshed a bar" % gmsh.option.getString("General.Version"))
            finally:
                gmsh.finalize()
        except Exception as e:
            failures.append("gmsh could not mesh anything from this interpreter (%s)" % e)

        if not failures:
            grid = subprocess.run(
                ["ElmerGrid", "14", "2", msh, "-autoclean", "-out", os.path.join(work, "mesh")],
                capture_output=True,
                text=True,
            )
            names = os.path.join(work, "mesh", "mesh.names")
            if grid.returncode != 0 or not os.path.isfile(names):
                failures.append("ElmerGrid could not convert a gmsh mesh: %s" % (grid.stdout + grid.stderr)[-500:])

        if not failures:
            ids = dict(re.findall(r"\$\s*(\w+)\s*=\s*(\d+)", open(names).read()))
            with open(os.path.join(work, "case.sif"), "w") as f:
                f.write("""
Header
  Mesh DB "." "mesh"
End
Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
End
Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(2) = 1 2
End
Solver 1
  Equation = Heat Equation
  Procedure = "HeatSolve" "HeatSolver"
  Variable = Temperature
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
End
Solver 2
  Exec Solver = After Simulation
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = "scalars.dat"
  Variable 1 = Temperature
  Operator 1 = max
End
Material 1
  Heat Conductivity = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = %s
  Temperature = 0
End
Boundary Condition 2
  Target Boundaries(1) = %s
  Temperature = 10
End
""" % (ids["cold"], ids["hot"]))
            solve = subprocess.run(["ElmerSolver", "case.sif"], cwd=work, capture_output=True, text=True)
            scalars = os.path.join(work, "scalars.dat")
            if solve.returncode != 0 or not os.path.isfile(scalars):
                failures.append("ElmerSolver could not solve a heat problem: %s" % (solve.stdout + solve.stderr)[-800:])
            else:
                hottest = float(open(scalars).read().split()[0])
                print("ElmerSolver solved a bar held at 0 and 10: hottest point %.3f" % hottest)
                if abs(hottest - 10.0) > 1e-6:
                    failures.append("ElmerSolver's answer to a bar held at 0 and 10 peaks at %g, not 10" % hottest)

    if failures:
        print("\nThis image cannot run the Elmer analyses:")
        for failure in failures:
            print("  * %s" % failure)
        return 1
    print("\nElmer and the mesher are both usable from this interpreter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
