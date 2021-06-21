```
      _   _ ______ __  __  ____             .                 
     | \ | |  ____|  \/  |/ __ \           ":"
     |  \| | |__  | \  / | |  | |        ___:____     |"\/"|
     | . ` |  __| | |\/| | |  | |      ,'        `.    \  / 
     | |\  | |____| |  | | |__| |      |  O        \___/  |
     |_| \_|______|_|  |_|\____/    -~^~^~^~^~^~^~^~^~^~^~^~^~
                                      -~^~^ ~^~ ^~^~^ ~^~^-

           ## pre-alpha ##           Mobilis adversus flumine. 
```

# nemo (pre-alpha version) #

Nemo is a lightweight, easy to understand (magneto)-hydrodynamics code
leveraging the robustness of 'Finite Volume' (FV) methods and the
efficiency of nodal 'Discontinuous Galerkin Spectral Element Methods' (DGSEM).

NOTE: The code is experimental and documentation is lacking.

## Installation ##

The following steps have to be done only once per machine.

```bash
# Clone nemo.
git clone https://github.com/jmark/nemo.git

# Remember where to find nemo.
export NEMO_HOME="$(pwd)/nemo"

# Change to the 'site' matching your machine.
cd "${NEMO_HOME}/sites/<my_machine>/<my_compiler>"

# Install external libraries. May take a little while.
bash install_deps.sh

# Remember where to find the external libraries.
export NEMO_INSTALL="$(pwd)/__install__'
```

You are now ready to roll! Happy coding!

## Quick Start ##

```bash
# Create and change to a working directory.
mkdir -p "/path/of/my/work/directory/"
cd "/path/of/my/work/directory/"

# Make a copy of the handy start-up script which matches your machine.
cp "${NEMO_HOME}/sites/<my_machine>/<my_compiler>/nemo.sh" .

# Open the script with your favorite text editor and follow the comments.
vim nemo.sh

# Build and run!
./nemo.sh all
```

If successful there are now a 'build' and an 'output' directory. Happy simulating!

## Third Parties ##

Nemo either uses or links to following third party softare:

    - fypp: https://github.com/aradi/fypp
    - HDF5: https://www.hdfgroup.org/solutions/hdf5/
    - p4est: https://www.p4est.org/
    - FMM3D: https://github.com/flatironinstitute/FMM3D

## Code Structure ##

The directory structure of the sources is simple and intuitively named. It contains:

```
source/
    - checkpoint/       ... I/O of checkpoints in the venerable FLASH file format.
    - equations/        ... Implementation of the various conservation laws.
    - kernel/           ... Solvers for the spatial part of the PDE.
    - mesh/             ... Mesh routines. Bindings to 'p4est'.
    - physics/          ... Source terms: (self-)gravity, (radiative) cooling, etc.
    - setups/           ... Definition of the simulation setups.
    - timedisc/         ... Time discretization, i.e. Runge-Kutta (RK) methods.
    - utils/            ... Auxiliary modules for general usage. Not specific to nemo.

    - nemo_prog.f90     ... Entry point to the code.
    - share_mod.f90     ... Global objects accessible by all modules, i.e. runtime, mesh, traceback.
    - driver_mod.f90    ... Game loop driving the simulation forward in time.
    - runtime_mod.f90   ... Definition of the 'runtime_t' type containing general
                            information such as the current simulation time, RK stage, etc.
```

All Fortran files implement **one** module and must follow the same naming
convention. That is the module name must end with '\_mod' and matches the file
name sans file extension (.f90). For example the file 'my\_module\_mod.f90' would look as
follows.
```fortran
module my_module_mod

!! ...

end module
```

Nemo offers a plethora of methods to choose from for meshing, hydrodynamics,
time integration and more. Different implementations of equal functionality
provide identical interfaces and are grouped under a common directory. There
are for example a variety of time integrators available which are arranged
in a hierarchy of sub-folders under 'source/timedisc/'.
```
source/timdisc/
    - euler/                ... first order Euler method
    - rungekutta/5-4        ... low storage 4th order RK
    - rungekutta/ssp/5-4    ... strong-stability preserving 4th order RK
    - ...
```

Ideally, one can exchange a method with another without any code modifications.

## Building Process ##

Nemo solely consists of Fortran (\*.f90) and Python (\*.py) files. There are no
configuration files read in at runtime. All configuration is done via Fortran
files and via compile-time pragmas (conditional compilation). That means any
changes in the setup definition are followed by an ordinary rebuild. This
is not a problem since building of nemo is slick and fast due to proper
module dependency resolution and parallelization.

The building of nemo consists of three phases: 'preprocessing', 'source code
filtering' and the actual 'compilation' with the Fortran compiler. The
following paragraphs briefly explain the three phases and show their respective
commands to evoke in a shell session on a typical Linux cluster.

### Phase 1: Preprocessing ###

Every simulation setup consists of at least one 'preproc.py' script. This
preprocessing script is a normal Python program which is supposed to be modified
by the user. The 'nemo' Python module provides a 'Makefile' object where one
registers all Fortran files needed for a successful compilation.

The Fortran files are parsed for compile-time pragmas and module dependencies.
In the end the script generates a standard-compliant Makefile which can be
called by Make.

This system enables the user to conveniently choose for example between
different hydrodynamics solvers or to activate source terms by mere registering
the respective Fortran modules. Think of assembling a recipe by enumerating all
needed ingredients.

The following code lists an exemplary snippet of such a 'preproc.py' script.
```python
#!/usr/bin/env python

import sys,os

home = os.environ['NEMO_HOME']
setup = os.environ['NEMO_SETUP']

sys.path.append(os.path.join(home,'scripts'))
from nemo import Makefile

mk = Makefile(home)

mk.add('..._mod.f90')

# ...

mk.generate()
```

Before using nemo it is advisable to change into a separate working directory
away from the nemo source tree. Hereby the nemo sources are not polluted with
build files and outputs from the simulation.

```bash
cd "/path/to/my/work/directory/"
```

The preprocessing script expects some environment variables in order to find
the nemo sources and to know which Fortran compiler to use.

```bash
export NEMO_HOME="/path/to/the/nemo/root/directory/"
export NEMO_SETUP="${NEMO_HOME}/source/setups/..."

export FC="gfortram" # mpif90, ifort, ...
export FCFLAGS="-O3"
```

The preprocessing is initiated with
```bash
python "${NEMO_SETUP}/preproc.py"
```

Per default the generated Makefile is written to a 'build/' directory which is
automatically created when it does not yet exist.

### Phase 2: Filtering ###

Nemo utilizes a sophisticated macro system provided by the source code filter
program 'fypp'. Fypp reads an original Fortan file from the nemo source tree
and writes a filtered version to the 'build/' directory.

The source code filtering is initiated with
```bash
make -j4 -C build/ filter
```

When finished a completely independent and specially tailored version of nemo
resides now in 'build/'.

### Phase 3: Compilation ###

The final phase is initiated by
```bash
make -j4 -C build/ compile
```

## Running ##
After successful compilation nemo is ready to roll!

```bash
./build/nemo
```
