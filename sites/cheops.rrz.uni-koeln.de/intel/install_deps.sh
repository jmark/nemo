#!/bin/bash -l

module purge
module load intel
module load intelmpi
module load git

export F9X=mpiifort
export CC=mpiicc
export CXX=mpicxx
export FC=mpiifort

set -eu -o pipefail

make -C __install__/p4wrap/sites/intel/

exit

mkdir -p __install__
pushd __install__

# P4WRAP
if test ! -e p4wrap/.done
then
    git clone 'https://github.com/jmark/p4wrap.git'
    pushd p4wrap
        git submodule update --init --recursive
        make -C sites/intel
        touch .done
    popd
fi

# FMM3D
if test ! -e FMM3D/.done
then
    git clone 'https://github.com/jmark/FMM3D.git'
    pushd FMM3D
        ln -s make.inc.icc make.inc
        make -j4 OMP=ON lib
        touch .done
    popd
fi

popd

echo
echo "!! SUCCESS !!"
