#!/usr/bin/env python

import sys,os

home = os.environ['NEMO_HOME']
setup = os.environ['NEMO_SETUP']
sys.path.append(os.path.join(home,'scripts'))

from nemo import Makefile

flags = dict(
    PP_N_NODES              = 8,
    PP_KERNEL_DGFV_QMAX     = 3,
    PP_MESH_PERIODIC        = 0,
    DO_PROFILING            = 1,
)

mk = Makefile(home,flags)

kernel = home + '/source/kernel/dgfv-smoothing/3d'
#kernel = home + '/source/kernel/fv-minmod/3d'

#timedisc    = home + '/source/timedisc/rungekutta/5-4'
timedisc    = home + '/source/timedisc/rungekutta/ssp/5-4'

equations = home + '/source/equations/euler/polytropic/3d'
#riemann = equations + '/riemann/rusanov'
#riemann = equations + '/riemann/hllc'

# mk.add(equations + '/riemann/roe-entropy-fix-no-carbuncle/riemann_mod.f90',alias='riemann_inner_mod')
# mk.add(equations + '/riemann/roe-entropy-fix-no-carbuncle/riemann_mod.f90',alias='riemann_boundary_mod')

mk.add(equations + '/riemann/hlle/riemann_mod.f90',alias='riemann_inner_mod')

mk.add(equations + '/riemann/eulerian-slip-wall/riemann_mod.f90',alias='riemann_north_mod')
mk.add(equations + '/riemann/eulerian-slip-wall/riemann_mod.f90',alias='riemann_west_mod')
mk.add(equations + '/riemann/eulerian-slip-wall/riemann_mod.f90',alias='riemann_front_mod')

mk.add(equations + '/riemann/hlle/riemann_mod.f90',alias='riemann_south_mod')
mk.add(equations + '/riemann/hlle/riemann_mod.f90',alias='riemann_east_mod')
mk.add(equations + '/riemann/hlle/riemann_mod.f90',alias='riemann_back_mod')

# mk.add(equations + '/riemann/rusanov/riemann_mod.f90',alias='riemann_inner_mod')
#mk.add(equations + '/riemann/hlle/riemann_mod.f90',alias='riemann_inner_mod')
#mk.add(equations + '/riemann/hlle+roe-entropy-fix/riemann_mod.f90',alias='riemann_boundary_mod')

#mk.add(equations + '/two_point/standard/two_point_flux_mod.f90')
#mk.add(equations + '/two_point/chandrashekar/two_point_flux_mod.f90')

mk.add(home + '/source/share_mod.f90')
mk.add(home + '/source/driver_mod.f90')
mk.add(home + '/source/runtime_mod.f90')

mk.add(home + '/source/utils/*.f90')
mk.add(home + '/source/mesh/**.f90')

mk.add(home + '/source/checkpoint/hdf5/flash/*.f90')

mk.add(equations + '/*.f90')
mk.add(kernel + '/*.f90')
mk.add(timedisc + '/*.f90')

mk.add(setup + '/../_common_/*_mod.f90')
mk.add(setup + '/*_mod.f90')

mk.generate()
