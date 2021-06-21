#!/usr/bin/env python

import sys,os

home = os.environ['NEMO_HOME']

sys.path.append(os.path.join(home,'scripts'))

from nemo import Makefile

mk = Makefile(
    root = home,
    #flags = dict(PP_N_NODES = 8, PP_KERNEL_DGFV_QMAX = 2, PP_SPLIT_FORM = 1, PP_MESH_PERIODIC = 0, DO_PROFILING = 1),
)

sourcedir   = home + '/source'
SETUPDIR    = sourcedir + '/setups/euler/3d/evrard'

hydro       = home + '/source/flux/euler/dgfv/3d'
#hydro       = home + '/source/flux/euler/dgfv/3d/multilevel-lobatto'
#hydro       = home + '/source/flux/euler/dgfv/3d/multilevel'
#hydro       = home + '/source/flux/euler/dgfv/3d/singlelevel'
#hydro       = home + '/source/flux/euler/dgfv/3d/singlelevel-multiflux'
kernel       = home + '/source/kernel/dgfv/3d'

#timedisc    = home + '/source/timedisc/rungekutta/5-4'
timedisc    = home + '/source/timedisc/rungekutta/ssp/5-4'

equations = home + '/source/equations/euler/polytropic/3d'
#riemann = equations + '/riemann/rusanov'
#riemann = equations + '/riemann/hllc'
#riemann = equations + '/riemann/hlle'

#mk.add(equations + '/two_point/standard/two_point_flux_mod.f90')
mk.add(equations + '/two_point/chandrashekar/two_point_flux_mod.f90')

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

mk.add(home + '/source/share_mod.f90')
mk.add(home + '/source/driver_mod.f90')
mk.add(home + '/source/runtime_mod.f90')

mk.add(home + '/source/common/*.f90')
mk.add(home + '/source/mesh/**.f90')

mk.add(home + '/source/checkpoint/base/*.f90')
mk.add(home + '/source/checkpoint/hdf5/p4est/*.f90')

#mk.add(riemann + '/*.f90')

mk.add(equations + '/riemann/riemann_hlle_mod.f90')
mk.add(equations + '/riemann/riemann_roe_entropy_fix_mod.f90')

mk.add(equations + '/*.f90')
mk.add(kernel,'*.f90')

mk.add(timedisc,'*.f90')

mk.add(SETUPDIR + '/../_common_/*_mod.f90')
mk.add(SETUPDIR + '/*_mod.f90')

mk.add(home + '/source/physics/self_gravity/*_mod.f90')

mk.generate()
