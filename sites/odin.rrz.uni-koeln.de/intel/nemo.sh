#!/bin/bash -l

# =========================================================================== #
# =========================================================================== #

# Where to find the root directory of nemo?
export NEMO_HOME="${HOME}/nemo"

# Where are the external libraries installed? (Probably already correct.)
export NEMO_INSTALL="${NEMO_HOME}/sites/odin.rrz.uni-koeln.de/intel/__install__"

# Which setup should nemo run?
export NEMO_SETUP="${NEMO_HOME}/source/setups/mhd/3d/cloud-shock"

# Which job name should slurm display?
export NEMO_JOBNAME="nemo: 3D Cloud-Shock"

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

    all         : execute above commands in order: preproc > filter > compile

    run         : run 'nemo' locally
    sbatch      : submit 'nemo' via sbatch

    all+run     : execute above commands in order: preproc > filter > compile > run
    all+sbatch  : execute above commands in order: preproc > filter > compile > sbatch
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

FCFLAGS += -lp4est
FCFLAGS += -lhdf5_fortran

# 'p4wrap' is stattically linked in.
FCFLAGS += -L$(NEMO_INSTALL)/p4wrap/sites/odin.rrz.uni-koeln.de/intel/build
FCFLAGS += -lp4wrap

# 'fmm3d' is stattically linked in.
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
module load p4est/2.2
module load hdf5/1.12.0

make -j6 -C build compile

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'run'
then

module purge
module load p4est/2.2
module load hdf5/1.12.0

rm -rf output
mkdir -p output

export OMP_NUM_THREADS="6"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_SCHEDULE=dynamic

pushd output > /dev/null
../build/nemo | tee out.log
popd > /dev/null

# =========================================================================== #
# =========================================================================== #
# Provide your credentials and set partition, number of nodes, etc.

elif test "${COMMAND}" == 'sbatch'
then

# Remove this line.
echo
echo "Provide your credentials in order to use sbatch." 1>&2; exit 1;

PARTITION="gassner"

NODES='8'
CPUS_PER_TASK='16'
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

#SBATCH --mem=30gb
#SBATCH --nodes=${NODES}
#SBATCH --ntasks=${NTASKS}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}

#SBATCH --time=24:00:00

#SBATCH --parsable
#SBATCH --error=out.log
#SBATCH --output=out.log

#SBATCH --job-name="${NEMO_JOBNAME}"

module purge
module load p4est/2.2
module load hdf5/1.12.0

export OMP_NUM_THREADS=${CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_SCHEDULE=dynamic

srun -n "\${SLURM_NTASKS}" ../build/nemo
EOT
)

rm -rf output
mkdir output

pushd 'output' > /dev/null
echo "${SBATCH_SCRIPT}" | sbatch | tee JOBID
popd > /dev/null

# =========================================================================== #
# =========================================================================== #
# Nothing to change here.

elif test "${COMMAND}" == 'all'
then

$0 preproc && $0 filter && $0 compile

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
