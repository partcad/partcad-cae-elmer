# The runtime the Elmer analyses run best in.
#
# It carries what pip cannot: Elmer, which is a native solver, and a gmsh for
# 64-bit ARM Linux, for which pip has no wheel. Everything else the scripts
# import is declared in `partcad.yaml` as `pythonRequirements` and installed by
# PartCAD into the environment it creates -- this image is a *runtime*, and the
# package has to work in a conda or venv sandbox on a host with Elmer installed
# too.
#
# It does NOT carry this package's code. PartCAD mounts the package at run time,
# so a user who fetches a newer version gets that version.
ARG PYTHON_VERSION=3.11
ARG ARCH=amd64

# Elmer is in neither Debian nor conda-forge, so it is built here, from a pinned
# release, in a stage of its own: the compilers and the source tree stay behind
# and what is copied out is the installed prefix, a few tens of megabytes.
#
# Without MPI and without the GUI -- one process per analysis, driven by a file
# -- and with OpenMP, which is how ElmerSolver uses more than one core in a
# single process. UMFPACK and ARPACK come with Elmer's own sources, and are what
# the direct solves and the modal analysis use.
FROM ghcr.io/partcad/partcad-container-python:py${PYTHON_VERSION}-${ARCH} AS build
USER root
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install --yes --no-install-recommends \
       git ca-certificates cmake make gfortran g++ libopenblas-dev liblapack-dev \
  && rm -rf /var/lib/apt/lists/*
ARG ELMER_RELEASE=release-26.2.1
RUN git clone --depth 1 --branch "${ELMER_RELEASE}" https://github.com/ElmerCSC/elmerfem.git /src/elmerfem \
  && cmake -S /src/elmerfem -B /src/build \
       -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/elmer \
       -DWITH_MPI=OFF -DWITH_OpenMP=ON -DWITH_ELMERGUI=OFF -DWITH_LUA=OFF \
       -DWITH_Mumps=OFF -DWITH_Hypre=OFF \
  && cmake --build /src/build -j"$(nproc)" \
  && cmake --install /src/build \
  && rm -rf /src

FROM ghcr.io/partcad/partcad-container-python:py${PYTHON_VERSION}-${ARCH}

LABEL partcad.image="1"
LABEL org.opencontainers.image.source="https://github.com/partcad/partcad-cae-elmer"

USER root
ENV DEBIAN_FRONTEND=noninteractive

# What the solver links against at run time, and the mesher. Debian builds
# `python3-gmsh` for amd64 and arm64 alike, which is the reason this image can
# promise gmsh on ARM at all.
RUN apt-get update \
  && { apt-get install --yes --no-install-recommends \
         libgfortran5 libgomp1 libopenblas0 liblapack3 python3-gmsh \
       || { echo "=== apt cannot install what this image is for ==="; \
            . /etc/os-release && echo "base: ${PRETTY_NAME}"; \
            apt-cache policy libgfortran5 libopenblas0 python3-gmsh || true; \
            exit 1; }; } \
  && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/elmer /opt/elmer

# Where every way of looking finds it: on PATH through /usr/local/bin, which
# `elmer_common.find_elmer()` searches even when PATH says nothing, and
# ELMER_HOME, which is how ElmerSolver finds its own solver modules.
RUN ln -s /opt/elmer/bin/ElmerSolver /usr/local/bin/ElmerSolver \
  && ln -s /opt/elmer/bin/ElmerGrid /usr/local/bin/ElmerGrid
ENV ELMER_HOME=/opt/elmer

# gmsh's Python API is one ctypes module wrapping `libgmsh.so`, installed by
# Debian for Debian's interpreter; the one in this image is python.org's. A
# ctypes wrapper works from either, so the one file is copied where both can see
# it and named on PYTHONPATH -- PYTHONPATH rather than site-packages, because
# PartCAD runs the scripts in a virtual environment it creates at run time.
RUN mkdir -p /opt/pc-site \
  && cp /usr/lib/python3/dist-packages/gmsh.py /opt/pc-site/gmsh.py
ENV PYTHONPATH=/opt/pc-site

# What a caller may run, extending the base image's allowlist rather than
# replacing it.
ENV PC_CONTAINER_ALLOWED_COMMANDS='{"python": "/usr/local/bin/python3", "ElmerSolver": "/usr/local/bin/ElmerSolver", "ElmerGrid": "/usr/local/bin/ElmerGrid"}'

# Two checks: the base image's own, for the contract PartCAD talks to a
# container through; and this package's, which solves a small heat problem with
# the solver and the mesher from *this* interpreter, the way the analyses will.
COPY verify_image.py /tmp/verify_image.py
RUN python3 /pc/verify.py \
  && python3 /tmp/verify_image.py \
  && rm /tmp/verify_image.py

USER pc
