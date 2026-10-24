# partcad-cae-elmer

Static stress, modal, steady heat transfer and laminar flow analysis for
[PartCAD](https://github.com/partcad/partcad), by [Elmer](https://www.elmerfem.org/). It implements all four of
PartCAD's analyses:

```shell
pc cae fea :bracket        # does it hold up under what it carries
pc cae dynamic :bracket    # at which frequencies does it ring, held as it is
pc cae thermal :heatsink   # how hot does it get, heated and cooled as it is
pc cae cfd :duct           # how does a fluid flow through it
pc test                    # every analysis a part declares, as a check that fails on any finding
```

Name it for one run with `-i //pub/feature/cae/elmer:<analysis>`, for one of a part's scenarios with
`implementation:` in it, or for everything with the `cae<Analysis>Implementation` user configuration option
(`caeDynamicImplementation` and `caeThermalImplementation` default to this package).

## What each analysis reads

Each analysis reads its own section of the part, a map of named scenarios: `pc cae <analysis>` runs each, and
`-f <scenario>` one. Every key names interfaces the part implements, and is applied to the **faces each port
covers**: PartCAD imprints a port's sketch on the part before this package runs, so an inlet drawn as a 3 mm
circle is a 3 mm circle; a port with no sketch covers the face it lies on.

| analysis | section keys | what they mean here |
| --- | --- | --- |
| `fea` | `fix:`, `load:`, `gravity:` | held still in every direction; a force spread over the faces, along the direction the part states (`[x, y, z]` in its frame) or along -Z for a magnitude; the acceleration of free fall the part is weighed under, in its frame |
| `dynamic` | `fix:` | held still; the lowest `modes` natural frequencies and their shapes are found |
| `thermal` | `temperature:`, `heat:`, `convection:` | held at a temperature (C); a power put in (W); cooled by a film coefficient (W/(m^2*K)) against `ambient_temperature`. Every other surface is insulated |
| `cfd` | `load:`, `outlet:`, `fix:` | an inlet pushed on with its force over its area; an outlet at the reference pressure; walls. Every surface that is neither an inlet nor an outlet is a no-slip wall |

```yaml
parts:
  cantilever:
    type: build123d
    implements:
      beam-root: [[-50, 0, 0], [0, 0, 1], 0]
      beam-tip: [[50, 0, 0], [0, 0, 1], 0]
    dynamic:
      clamped:
        fix:
          - beam-root
    thermal:
      heated:
        temperature:
          beam-root: 20 C
        heat:
          beam-tip: 1 W
```

### What the numbers mean

A `load` value is a force, and PartCAD turns it into newtons before this package sees it. Write a number and a
unit -- `n`, `mn`, `kn` or `newton` for a force, `mg`, `g`, `kg`, `ton`, `tonne`, `lb` or `pound` for a mass,
matched case-insensitively, with or without a space and a plural `s` -- or a **bare number, which is a mass in
kilograms**, weighed at 9.8 N/kg: `hook: 5` is 49 N. A direction is three of those, each read the same way, so
`[0, 0, -5]` is the same 49 N straight down, and `[0 N, -2.7 N, -155.2 N]` is newtons as written. What a
simulation measured is in newtons, so write the unit: `[0, -2.7, -155.2]` would be read as 1.5 kN.
`cfd` reads its `load:` the same way: the newtons that push the fluid in, spread over the inlet's area.

`gravity:` is an acceleration, m/s² along each of the part's axes -- a number, or a string saying so
(`9.81 m/s^2`). The part's own weight is its mass, `density` times its volume, under that acceleration, in
place of the 9.81 m/s² along -Z that `self_weight` applies otherwise; it is never a `load`.

What comes back is a `<part>.<analysis>.<scenario>.glb` the Digital Thread shows in the analysis's tab --
coloured by von Mises stress, by how far each point moves in the first mode, by temperature, or by pressure --
a measurement line in the log, and **findings**, which fail `pc test`:

* `fea`: a stress past the yield strength (an error) or past it divided by `safety_factor`, a displacement
  past `max_displacement`;
* `dynamic`: a first natural frequency under `min_frequency`, when that is set;
* `thermal`: a temperature past `max_temperature`, when that is set;
* `cfd`: a flow that did not settle, a Reynolds number past 2300 (where this laminar analysis stops
  applying), a speed past `max_velocity`, a flow rate short of `min_flow_rate` when that is set;
* every analysis: a condition whose port covers no surface of the part, and so did nothing.

### The pipeline

gmsh meshes the part into tetrahedra -- quadratic for the structural and thermal analyses, linear for the
flow -- with each condition's faces a named group (PartCAD's `cae_support.mesh()`). ElmerGrid converts the
mesh and says which number each group became. A solver input file is written for StressSolve (static, or as
an eigen analysis), HeatSolve or FlowSolve, and ElmerSolver runs it. The fields come back as an ASCII VTU
file, and the eigenvalues from the log.

### Measured

`examples/feature_cae` in the PartCAD repository has a cantilever and a pipe whose answers are formulas.
Through `pc cae`, in the image:

| analysis | closed form | this package |
| --- | --- | --- |
| `fea`: tip deflection, 100 N on a 100x10x10 mm steel bar | 0.2016 mm (with shear) | 0.2005 mm |
| `dynamic`: first bending frequency | 815.4 Hz (Euler-Bernoulli) | 813.25 Hz, twice -- the square section's two bending planes |
| `thermal`: hottest point, 1 W at the tip, root at 20 C | 40.0 C | 40.0 C |
| `cfd`: flow rate, 5 mN through a 10x100 mm bore of oil | 1.5625 mL/s | 1.551 mL/s |

Code_Aster, on the same mesh, gives the same displacement, frequencies and temperatures to every printed digit.

## Parameters

Set in a `cae: <analysis>:` section of the package that uses this one, or of the part, or in one scenario's
`params:`. See `partcad.yaml` for every default and why it is what it is.

| parameter | analyses | meaning |
| --- | --- | --- |
| `mesh_size` | all | largest tetrahedron, as a fraction of the part's largest dimension |
| `mesh_order` | fea, dynamic, thermal | 2 for quadratic tetrahedra, 1 for linear |
| `youngs_modulus`, `poissons_ratio`, `density` | fea, dynamic | the material, in Pa and kg/m^3 |
| `yield_strength`, `safety_factor`, `max_displacement`, `self_weight` | fea | the thresholds, and gravity (9.81 m/s² along -Z, unless the part states `gravity:`) |
| `modes`, `min_frequency` | dynamic | how many modes, and the lowest first frequency allowed (Hz) |
| `conductivity`, `ambient_temperature`, `max_temperature` | thermal | W/(m*K), and Celsius |
| `density`, `viscosity`, `iterations`, `max_velocity`, `min_flow_rate` | cfd | the fluid (kg/m^3, Pa*s), the solve, the thresholds |

## What it needs

**A container runtime, or Elmer installed.**

With a container runtime PartCAD runs the analyses in the image `partcad.yaml` names: PartCAD's Python base
image, with Elmer built from a pinned release and Debian's `python3-gmsh`, for amd64 and arm64. Nothing to
install.

Without one, PartCAD runs them in a conda or venv sandbox on the host, and Elmer has to be there.
`ElmerSolver` and `ElmerGrid` are found, in this order, in the directory `PARTCAD_ELMER_BIN` names, in
`$ELMER_HOME/bin`, on `PATH`, and in `/opt/elmer/bin`, `/usr/local/bin`, `/usr/bin`, `/opt/local/bin`.
ElmerSolver finds its solver modules through `ELMER_HOME`, which is derived from where the executable is when
it is not set. Ways to have it: CSC's Ubuntu packages (`ppa:elmer-csc-ubuntu/elmer-csc-ppa`, package
`elmerfem-csc`), the Windows installer from elmerfem.org, or a build of
[ElmerCSC/elmerfem](https://github.com/ElmerCSC/elmerfem) -- the Dockerfile is a working recipe.

The daemon is what runs the analysis, so a variable set in your shell reaches it only when the daemon is
started from that shell: `pc daemon stop` first.

gmsh comes from pip everywhere pip has a wheel, which is everywhere except 64-bit ARM Linux; there, use the
image.

### The pin

`dockerImage:` names the image by a hash of what it is built from -- the Dockerfile, `verify_image.py` and the
Python version -- not by a version. Run `./image-tag.sh` after changing any of them and update every pin; the
`Image` workflow refuses to build if they disagree, and publishes `<tag>-amd64` and `<tag>-arm64` on a push to
`main`. The arm64 build compiles Elmer under emulation and takes a while.

## Tests

```shell
pytest
```

beside an installed `partcad` (the scripts import PartCAD's sandbox helpers; `conftest.py` finds them). The
readers, the input files and the messages a machine without Elmer gets are checked everywhere; the example
cantilever is analysed three ways and compared with its closed forms where Elmer, gmsh and build123d are all
present.

To keep a case for inspection, set `PARTCAD_ELMER_CASE` to an empty directory before the daemon starts.

## License

Apache License 2.0. Elmer's solver library is LGPL-2.1 and ElmerGrid GPL-2.0; neither is part of this
package. The image builds them from CSC's published sources at the release the Dockerfile names.
