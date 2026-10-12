#
# PartCAD, 2026
#
# Licensed under Apache License, Version 2.0.
#
"""Finding Elmer, running it, and reading what it wrote.

Everything the four analysis scripts share that is about Elmer rather than
about the physics: where `ElmerSolver` and `ElmerGrid` are, how a gmsh mesh
becomes an Elmer mesh directory, how a solver input file (`.sif`) is run, and
how its results come back -- the scalars it was asked to save, the eigenvalues
it printed, and the fields it wrote as an ASCII VTU file.

**Where it is.** `ElmerSolver` and `ElmerGrid` are looked for in the directory
`PARTCAD_ELMER_BIN` names, then in `$ELMER_HOME/bin`, then on `PATH`, then where
the usual installations put them. ElmerSolver finds its own solver modules
through `ELMER_HOME`, so when that is not set it is derived from where the
executable was found -- `<prefix>/bin/ElmerSolver` has its modules under
`<prefix>/share/elmersolver`.
"""

import os
import platform
import re
import shutil
import subprocess

# The meta-wrapper's directory is on 'sys.path' when PartCAD runs a script, and
# 'cae_support' lives there: what every implementation does around its solver.
import cae_support as support  # noqa: E402

BIN_ENV = "PARTCAD_ELMER_BIN"
DIRECTORIES = (
    "/opt/elmer/bin",
    "/usr/local/bin",
    "/usr/bin",
    "/opt/local/bin",
    "C:\\Program Files\\Elmer 9.0-Release\\bin",
)

_found = None


def find_elmer():
    """The paths of ElmerSolver and ElmerGrid, and the environment to run them in."""
    global _found
    if _found is not None:
        return _found

    def executable(directory, name):
        for candidate in (name, name + ".exe"):
            path = os.path.join(directory, candidate)
            if os.path.isfile(path) and os.access(path, os.X_OK):
                return path
        return None

    places = []
    named = os.environ.get(BIN_ENV)
    if named:
        if not os.path.isdir(named):
            raise support.ToolMissing("%s points at something that is not a directory: %s" % (BIN_ENV, named))
        places.append(named)
    if os.environ.get("ELMER_HOME"):
        places.append(os.path.join(os.environ["ELMER_HOME"], "bin"))
    on_path = shutil.which("ElmerSolver")
    if on_path:
        places.append(os.path.dirname(os.path.realpath(on_path)))
    places.extend(DIRECTORIES)

    for directory in places:
        solver = executable(directory, "ElmerSolver")
        grid = executable(directory, "ElmerGrid")
        if solver and grid:
            env = dict(os.environ)
            if not env.get("ELMER_HOME"):
                env["ELMER_HOME"] = os.path.dirname(os.path.dirname(os.path.realpath(solver)))
            _found = (solver, grid, env)
            return _found

    raise support.ToolMissing(
        "Elmer is not available on this %s-%s machine. Two ways to have it: run with a container runtime "
        "available, which gets the image named by 'dockerImage:' in partcad.yaml and carries Elmer -- or "
        "install it here: CSC's packages ('add-apt-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa' and "
        "'apt install elmerfem-csc' on Ubuntu), the Windows installer from elmerfem.org, or a build of "
        "https://github.com/ElmerCSC/elmerfem; then put ElmerSolver and ElmerGrid on PATH, set ELMER_HOME, "
        "or point %s at the directory holding them. Searched: %s"
        % (platform.system(), platform.machine(), BIN_ENV, ", ".join(places))
    )


def run(work, executable, *arguments, log):
    """Run one Elmer program in 'work' and return what it printed, kept as 'log'."""
    _solver, _grid, env = find_elmer()
    path = os.path.join(work, log)
    with open(path, "w") as f:
        completed = subprocess.run([executable, *arguments], cwd=work, stdout=f, stderr=subprocess.STDOUT, env=env)
    text = open(path, errors="replace").read()
    if completed.returncode != 0 or "ERROR::" in text:
        errors = [line for line in text.splitlines() if "ERROR" in line]
        tail = "\n".join((errors or text.strip().splitlines())[-15:])
        raise support.ToolFailed(
            "%s failed (exit %d):\n%s" % (os.path.basename(executable), completed.returncode, tail)
        )
    return text


def convert(work, msh):
    """Turn a gmsh mesh into an Elmer mesh directory; return its boundary numbers by name.

    ElmerGrid keeps gmsh's physical groups, numbering them afresh and writing
    the names it was given to 'mesh.names' -- which is the only way to know which
    number a group became.
    """
    _solver, grid, _env = find_elmer()
    run(work, grid, "14", "2", msh, "-autoclean", "-out", "mesh", log="log.ElmerGrid")
    names = os.path.join(work, "mesh", "mesh.names")
    if not os.path.isfile(names):
        raise support.ToolFailed("ElmerGrid wrote no mesh.names, so no boundary can be found by name")
    text = open(names).read()
    bodies, boundaries = {}, {}
    target = bodies
    for line in text.splitlines():
        if "boundaries" in line:
            target = boundaries
        match = re.match(r"\$\s*(\S+)\s*=\s*(\d+)", line)
        if match:
            target[match.group(1)] = int(match.group(2))
    return bodies, boundaries


def solve(work, sif):
    """Write a solver input file into 'work', run ElmerSolver on it, return the log."""
    solver, _grid, _env = find_elmer()
    with open(os.path.join(work, "case.sif"), "w") as f:
        f.write(sif)
    return run(work, solver, "case.sif", log="log.ElmerSolver")


def scalars(work, name="scalars.dat"):
    """The values a SaveScalars solver wrote, in the order it was asked for them."""
    path = os.path.join(work, name)
    if not os.path.isfile(path):
        raise support.ToolFailed("ElmerSolver saved no scalars: the solve did not finish")
    return [float(value) for value in open(path).read().split()]


def eigenvalues(log):
    """The eigenvalues an eigen analysis printed, smallest first.

    Elmer prints them as 'EigenSolve: <n>: <real> <imaginary>' under a line
    saying how many it computed; for an undamped structure they are real and
    are the squares of the angular frequencies.
    """
    values = []
    for match in re.finditer(r"EigenSolve:\s+(\d+):\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)", log):
        values.append(float(match.group(2)))
    return values


SIF_LINEAR = """
  Linear System Solver = Direct
  Linear System Direct Method = Umfpack
"""

SIF_ITERATIVE = """
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStabL
  BiCGStabl Polynomial Degree = 4
  Linear System Preconditioning = ILU1
  Linear System Max Iterations = 5000
  Linear System Convergence Tolerance = 1.0e-10
"""


def result_output(index):
    """Solver number 'index': one that writes every field, at every node, as an ASCII VTU file.

    ASCII so that it can be read back here with a regular expression rather than
    a VTK library; the bulk elements only, since the surface drawn is worked out
    from them.
    """
    return """
Solver %d
  Exec Solver = After Simulation
  Equation = ResultOutput
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = "result"
  Vtu Format = Logical True
  Ascii Output = Logical True
  Save Bulk Only = Logical True
  Save Geometry Ids = Logical False
End
""" % index


# --------------------------------------------------------------------------- #
# Reading a VTU file                                                          #
# --------------------------------------------------------------------------- #

# VTK cell type -> how many of a cell's points are corners, for the cells a 3D
# Elmer result is made of: linear and quadratic tetrahedra, hexahedra, prisms.
CORNERS = {10: 4, 24: 4, 12: 8, 25: 8, 13: 6, 26: 6, 14: 5}


def find_vtu(work):
    """The VTU file the result output wrote, wherever this Elmer version put it."""
    found = []
    for directory, _subdirectories, files in os.walk(work):
        found.extend(os.path.join(directory, name) for name in files if name.endswith(".vtu"))
    if not found:
        raise support.ToolFailed("ElmerSolver wrote no result file")
    return sorted(found)[-1]


def _arrays(section):
    """Name -> (components, values) for every ASCII DataArray in a piece of VTU."""
    arrays = {}
    for attributes, body in re.findall(r"<DataArray([^>]*)>(.*?)</DataArray>", section, re.S):
        # Case-insensitively: Elmer spells it 'FORMAT=' on one array of the
        # very file it spells 'format=' everywhere else.
        if 'format="ascii"' not in attributes.lower():
            raise support.ToolFailed("ElmerSolver wrote a binary VTU file; this reads only ASCII ones")
        name = re.search(r'Name="([^"]*)"', attributes)
        components = re.search(r'NumberOfComponents="(\d+)"', attributes)
        arrays[name.group(1) if name else ""] = (int(components.group(1)) if components else 1, body.split())
    return arrays


def _section(text, tag):
    match = re.search(r"<%s[^>]*>(.*?)</%s>" % (tag, tag), text, re.S)
    return match.group(1) if match else ""


def read_vtu(path):
    """An ASCII VTU file as (points, cells by type, point fields).

    'points' is a list of (x, y, z); 'cells' maps a VTK type to the point index
    lists of the cells of that type; 'fields' maps an array name to its values,
    one tuple per point.
    """
    text = open(path, errors="replace").read()
    ((_size, raw),) = _arrays(_section(text, "Points")).values()
    values = [float(value) for value in raw]
    points = [tuple(values[k : k + 3]) for k in range(0, len(values), 3)]

    topology = _arrays(_section(text, "Cells"))
    connectivity = [int(value) for value in topology["connectivity"][1]]
    offsets = [int(value) for value in topology["offsets"][1]]
    types = [int(value) for value in topology["types"][1]]
    cells = {}
    start = 0
    for end, kind in zip(offsets, types):
        cells.setdefault(kind, []).append(connectivity[start:end])
        start = end

    fields = {}
    for name, (size, raw) in _arrays(_section(text, "PointData")).items():
        numbers = [float(value) for value in raw]
        fields[name] = [tuple(numbers[k : k + size]) for k in range(0, len(numbers), size)]
    return points, cells, fields


def surface(points, cells):
    """The outer triangles of the volume cells in a VTU result."""
    triangles = []
    for kind, members in cells.items():
        if kind in CORNERS:
            triangles.extend(support.outer_triangles(members, CORNERS[kind]))
    return triangles


def magnitude(values):
    """The length of each vector in a field."""
    return [sum(component * component for component in value) ** 0.5 for value in values]


# --------------------------------------------------------------------------- #
# What every analysis says the same way                                       #
# --------------------------------------------------------------------------- #

ERROR = "error"
WARNING = "warning"
INFO = "info"

M2_PER_MM2 = 1e-6


def nothing_covered(record):
    """The finding for a condition whose port covers no surface of the part."""
    return {
        "severity": WARNING,
        "message": (
            "The port covers no surface of the part, so the condition on it did nothing. Place the port on "
            "the face it means, or draw its interface with a sketch that lies on one."
        ),
        "where": record.get("port") or record.get("interface_label") or "a port",
    }


def mesh_part(request, groups, work, order, rest=None):
    """Mesh the part with gmsh and convert it; return the gmsh mesh and Elmer's numbering.

    The boundary numbers ElmerGrid gave each named group are what a solver input
    file targets, so they come back beside the mesh: '(mesh, boundary ids by name)'.
    """
    msh = os.path.join(work, "part.msh")
    mesh = support.mesh(
        request, groups, size=float(request.get("mesh_size", 0.05)), order=order, rest=rest, files=[msh]
    )
    _bodies, boundaries = convert(work, msh)
    return mesh, boundaries


def targets(boundaries, names):
    """'Target Boundaries(n) = ...' for the groups 'names', as ElmerGrid numbered them."""
    ids = [boundaries[name] for name in names if name in boundaries]
    return "Target Boundaries(%d) = %s" % (len(ids), " ".join(str(i) for i in ids))


def header():
    """The opening of every solver input file: the mesh, the simulation, the one body."""
    return """
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
  Body Force = 1
End
"""
