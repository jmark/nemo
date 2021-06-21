#!/usr/bin/env python

import sys,os

root = os.environ['NEMO_ROOT']
site = os.environ['NEMO_SITE']

sys.path.append(os.path.join(root,'scripts'))

from nemo import Makefile

mk = Makefile(
    root = root,
    sites = site,
    flags = dict(PP_N_NODES = 4, PP_KERNEL_DGFV_QMAX = 2, PP_SPLIT_FORM = 0, PP_MESH_PERIODIC = 0),
)

sourcedir   = root + '/source'
SETUPDIR    = sourcedir + '/setups/euler/3d/self-gravity'

hydro       = root + '/source/flux/euler/dgfv/3d'
#hydro       = root + '/source/flux/euler/dgfv/3d/multilevel-lobatto'
#hydro       = root + '/source/flux/euler/dgfv/3d/multilevel'
#hydro       = root + '/source/flux/euler/dgfv/3d/singlelevel'
#hydro       = root + '/source/flux/euler/dgfv/3d/singlelevel-multiflux'
kernel       = root + '/source/kernel/dgfv/3d'

#timedisc    = root + '/source/timedisc/rungekutta/5-4'
timedisc    = root + '/source/timedisc/rungekutta/ssp/5-4'

equations = root + '/source/equations/euler/polytropic/3d'
riemann = equations + '/riemann/rusanov'
#riemann = equations + '/riemann/hllc'

#riemann     = hydro + '/riemann/standard/roe'
#riemann     = hydro + '/riemann/standard/roe-entropy-fix'
#riemann     = hydro + '/riemann/standard/hlle+roe-entropy-fix'
#riemann     = hydro + '/riemann/standard/hll'
#riemann     = hydro + '/riemann/standard/hllc'
#riemann     = hydro + '/riemann/standard/hllm'
#riemann     = hydro + '/riemann/standard/entropy-stable'

#riemann     = hydro + '/riemann/multidim/rusanov'
#riemann     = hydro + '/riemann/multidim/chandrashekar'

#twopoint    = hydro + '/equations/polytropic/two_point_flux/standard'
#twopoint    = hydro + '/equations/polytropic/two_point_flux/kennedy-gruber'
#twopoint    = hydro + '/equations/polytropic/two_point_flux/chandrashekar'
#twopoint    = hydro + '/equations/polytropic/two_point_flux/schnuecke'
#mk.add(twopoint + '/two_point_flux_mod.f90')

mk.add(root + '/source/share_mod.f90')
mk.add(root + '/source/driver_mod.f90')
mk.add(root + '/source/runtime_mod.f90')

mk.add(root + '/source/common/*.f90')
mk.add(root + '/source/mesh/*.f90')
mk.add(root + '/source/mesh/bindings/*.f90')

mk.add(root + '/source/checkpoint/base/*.f90')
mk.add(root + '/source/checkpoint/hdf5/p4est/*.f90')

mk.add(riemann + '/*.f90')
mk.add(equations + '/*.f90')
mk.add(kernel,'*.f90')

mk.add(timedisc,'*.f90')

mk.add(SETUPDIR + '/modules/*_mod.f90')

# mk.add(SETUPDIR + '/modules/setup_amr_mod.f90')
# mk.add(SETUPDIR + '/modules/setup_amr_sensor_mod.f90')
# mk.add(SETUPDIR + '/modules/setup_amr_types_mod.f90')

mk.generate()
