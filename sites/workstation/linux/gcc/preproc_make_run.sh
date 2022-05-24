#!/bin/bash -l

COMMAND="$1"; shift;

USAGE=$(cat << EOT
preproc & compile & run script for 'nemo'

USAGE:
    $(basename "$0") <COMMAND>

COMMANDS:
    preproc : generate 'Makefile'
    make    : run 'make'
    run     : start 'nemo'
    all     : execute above commands at once: preproc > make > run
EOT
)

# =========================================================================== #
# =========================================================================== #

# FIXME
echo "Modify script beforehand! Carefully read every line and adapt accordingly!" 1>&2; exit 1;

export NEMO_HOME="${HOME}/projects/nemo"
export NEMO_SITE="${NEMO_HOME}/sites/workstation/linux/gcc"
export NEMO_SETUP="${NEMO_HOME}/source/setups/euler/3d/boss-bodenheimer"

# =========================================================================== #
# =========================================================================== #

if test "${COMMAND}" == 'preproc'
then

FCFLAGS=$(cat << 'EOT'
FC  = mpif90

FCFLAGS =

#FCFLAGS += -g
FCFLAGS += -O3

FCFLAGS += -cpp
FCFLAGS += -pedantic
FCFLAGS += -pedantic-errors
FCFLAGS += -fall-intrinsics

FCFLAGS += -Wall
FCFLAGS += -Wno-unused-label
FCFLAGS += -Wno-unused-variable
FCFLAGS += -Wno-unused-function
FCFLAGS += -Wno-unused-dummy-argument
FCFLAGS += -Wno-maybe-uninitialized
FCFLAGS += -Wno-error=conversion
FCFLAGS += -Werror

FCFLAGS += -std=f2008
FCFLAGS += -ffree-form
FCFLAGS += -fimplicit-none
FCFLAGS += -ffree-line-length-none

#FCFLAGS += -fdefault-real-8
#FCFLAGS += -fcheck=array-temps

# enable this to compile with OpenMP
FCFLAGS += -fopenmp

FCFLAGS += -lhdf5_fortran
FCFLAGS += -lz

FCFLAGS += -L$(NEMO_SITE)/p4wrap/sites/gcc/build
FCFLAGS += -lp4wrap

FCFLAGS += -L$(NEMO_SITE)/FMM3D/lib-static
FCFLAGS += -lfmm3d
EOT
)

echo "${FCFLAGS}" | "${NEMO_SETUP}/preproc.py"

# =========================================================================== #
# =========================================================================== #

elif test "${COMMAND}" == 'make'
then

make -j -C build "$@"

# =========================================================================== #
# =========================================================================== #

elif test "${COMMAND}" == 'run'
then

CPU_CORES="$(grep -i 'cpu cores' /proc/cpuinfo | head -n1 | awk '{print $NF}')"
export OMP_NUM_THREADS="${CPU_CORES}"
#export OMP_NUM_THREADS="3"

#export OMP_PLACES=cores
#export OMP_PLACES=threads
#export OMP_PLACES=socket
#export OMP_PLACES='{0:3}'
#export OMP_PLACES='{2:3}'
export OMP_PLACES='{0:6}'

#export OMP_PROC_BIND=spread
#export OMP_PROC_BIND=true
export OMP_PROC_BIND=close

rm -rf output
mkdir -p output
pushd output
#mpirun -n 2 ../build/nemo
#mpirun ../build/nemo
../build/nemo
popd

# =========================================================================== #
# =========================================================================== #

elif test "${COMMAND}" == 'all'
then

$0 preproc && $0 make && $0 run

# =========================================================================== #
# =========================================================================== #

else

echo "${USAGE}" 1>&2; exit 1;

fi
