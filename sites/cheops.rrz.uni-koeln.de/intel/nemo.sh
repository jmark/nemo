#!/bin/bash -l

# =========================================================================== #
# =========================================================================== #

# Where to find the root directory of nemo?
export NEMO_HOME="${HOME}/projects/nemo"

# Where are the external libraries installed? (Probably already correct.)
export NEMO_INSTALL="${NEMO_HOME}/sites/cheops.rrz.uni-koeln.de/intel/__install__"

# Which setup should nemo run?
export NEMO_SETUP="${NEMO_HOME}/source/setups/euler/2d/sedov"

# Which job name should slurm display?
export NEMO_JOBNAME="nemo: 2D Sedov Blast"

# Remove this line.
echo "Modify the script first! Carefully read every comment and edit accordingly!" 1>&2; exit 1;

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

COMMAND="$1"; shift;
OUTPUTDIR="${1:-output}"; shift;

USAGE=$(cat << EOT
build & run script for 'nemo'

USAGE:
    $(basename "$0") <COMMAND>

COMMANDS:
    preproc     : generate 'Makefile'
    filter      : run source code filter
    compile     : compile to exectuable

    run         : run 'nemo' locally
    sbatch      : forward 'nemo' to sbatch

    all+run     : execute above commands in order: preproc > filter > compile > run
    all+sbatch  : execute above commands in order: preproc > filter > compile > sbatch

    all         : alias of 'all+run'
EOT
)

# =========================================================================== #
# =========================================================================== #
# Adapt compiler flags. Given defaults should already work just fine, though.

if test "${COMMAND}" == 'preproc'
then

FCFLAGS=$(cat << 'EOT'
# Makefile syntax

FC = mpiifort

FCFLAGS =

#FCFLAGS += -g
FCFLAGS += -O3

FCFLAGS += -r8 
FCFLAGS += -i4
FCFLAGS += -cpp
FCFLAGS += -fPIC
FCFLAGS += -real_size 64
FCFLAGS += -implicitnone
#FCFLAGS += -std08
FCFLAGS += -diag-disable 5462

# Enable OpenMP.
FCFLAGS += -qopenmp

FCFLAGS += -lz
FCFLAGS += -lm
FCFLAGS += -lblas
FCFLAGS += -llapack 
FCFLAGS += -lhdf5_fortran

FCFLAGS += -L$(NEMO_INSTALL)/p4wrap/sites/intel/build
FCFLAGS += -lp4wrap

FCFLAGS += -L$(NEMO_INSTALL)/FMM3D/lib-static
FCFLAGS += -lfmm3d
EOT
)

echo "${FCFLAGS}" | "${NEMO_SETUP}/preproc.py"

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'filter'
then

make -j6 -C build filter

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'compile'
then

module purge

module load intel/18.0
module load intelmpi/2018
module load hdf5/1.10.5-intelmpi

export OMP_NUM_THREADS=${CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export CC=mpiicc
export CXX=mpicxx
export FC=mpiifort

make -j6 -C build compile

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'run'
then

module purge

module load intel/18.0
module load intelmpi/2018
module load hdf5/1.10.5-intelmpi

export F9X=mpiifort
export CC=mpiicc
export CXX=mpicxx
export FC=mpiifort

rm -rf output
mkdir -p output/checkpoints

export OMP_NUM_THREADS="6"

pushd output
../build/nemo | tee out.log
popd

# =========================================================================== #
# =========================================================================== #
# Provide your credentials and set partition, number of nodes, etc.

elif test "${COMMAND}" == 'sbatch'
then

# Remove this line.
echo
echo "Provide your credentials in order to use sbatch." 1>&2; exit 1;

PARTITION="gassner-exclusive"

NODES='4'
CPUS_PER_TASK='24'
NTASKS_PER_NODE='1'
NTASKS="$(echo "${NODES} * ${NTASKS_PER_NODE}" | bc)"

cat << EOT

  setup:              ${SETUP}
  nodes:              ${NODES}
  cpus-per-task:      ${CPUS_PER_TASK}
  ntasks-per-node:    ${NTASKS_PER_NODE}
  
  ntasks:             ${NTASKS}

EOT

SBATCH_SCRIPT=$(cat << EOT
#!/bin/bash -l

#SBATCH --mail-type=all

#SBATCH --partition=${PARTITION}

#SBATCH --mem=40gb
#SBATCH --nodes=${NODES}
#SBATCH --ntasks=${NTASKS}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}

#SBATCH --time=192:00:00

#SBATCH --parsable
#SBATCH --error=out.log
#SBATCH --output=out.log

##SBATCH --workdir='output'
#SBATCH --job-name="${NEMO_JOBNAME}"

module purge

# Default crashes.
# module load intel
# module load intelmpi

# Older version works!
module load intel/18.0
module load intelmpi/2018
module load hdf5/1.10.5-intelmpi

export OMP_NUM_THREADS=${CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

srun -n "\${SLURM_NTASKS}" ../build/nemo
EOT
)

# export SLURM_ACCOUNT="AG..."
# export SLURM_MAILUSER="${USER}@uni-koeln.de"

# Alternatively put above two lines into a separate file
# and source the file in.
# source "${HOME}/slurm-user-config.sh"

rm -rf output
mkdir -p output/checkpoints

pushd 'output'
echo "${SBATCH_SCRIPT}" | sbatch \
    --account="${SLURM_ACCOUNT}" \
    --mail-user="${SLURM_MAILUSER}" | tee JOBID
popd

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'all'
then

$0 all+run

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'all+run'
then

$0 preproc && $0 filter && $0 compile && $0 run

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'all+sbatch'
then

$0 preproc && $0 filter && $0 compile && $0 sbatch

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

else

echo "${USAGE}" 1>&2; exit 1;

fi
