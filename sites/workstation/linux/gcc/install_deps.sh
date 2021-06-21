#!/bin/bash -l

# Further dependencies:
#   - gcc & gfortran
#   - MPI library: mpicc & mpif90
#   - hdf5 library with Fortran bindings

set -eu -o pipefail

mkdir -p __install__
pushd __install__

# P4WRAP
if test ! -e p4wrap/.done
then
    # git clone 'https://github.com/jmark/p4wrap.git'
    pushd p4wrap/sites/workstation/gcc
        # git submodule update --init --recursive
        # git clone https://github.com/cburstedde/p4est.git p4est
        # pushd p4est
        # git submodule init
        # git submodule update
        # popd
        make
        touch .done
    popd
fi

# FMM3D
if test ! -e FMM3D/.done
then
    # If you have problems compilling, try to add '-fallow-argument-mismatch'
    # to the compiler flags.

    git clone 'https://github.com/jmark/FMM3D.git'
    pushd FMM3D
        make -j OMP=ON lib
        touch .done
    popd
fi

# HDF5
if test ! -e hdf5-1.10.6/.done
then
    wget -c 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz'
    tar xf 'hdf5-1.10.6.tar.gz'

    pushd hdf5-1.10.6
        ./configure --enable-fortran --enable-parallel --prefix="$(pwd)/__prefix__"
        make -j8 && make install
        touch .done
    popd
fi

popd

echo
echo "!! SUCCESS !!"
