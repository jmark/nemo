#!/usr/bin/env python

import sys,os

home = os.environ['NEMO_HOME']
setup = os.environ['NEMO_SETUP']
sys.path.append(os.path.join(home,'scripts'))

from nemo import Makefile

flags = dict(
    PP_N_NODES          = 8,
    PP_KERNEL_DGFV_QMAX             = 3,
    PP_MESH_PERIODIC    = 1,
    DO_PROFILING        = 0,
    PP_SPLIT_FORM       = 0,
)

mk = Makefile(home,flags)

#kernel = home + '/source/kernel/dgfv/2d'
kernel = home + '/source/kernel/dgfv-smoothing/2d'
#kernel = home + '/source/kernel/fv-mc/2d'
#kernel = home + '/source/kernel/fv/2d'

#kernel = home + '/source/kernel/fv-minmod/2d'
#kernel = home + '/source/kernel/fv/2d'
# timedisc    = home + '/source/timedisc/rungekutta/5-4'

timedisc = home + '/source/timedisc/rungekutta/ssp/5-4'

equations = home + '/source/equations/mhd/polytropic/2d'

mk.add(equations + '/riemann/rusanov/riemann_mod.f90',alias='riemann_inner_mod')

for d in 'north south west east'.split():
    mk.add(equations + '/riemann/rusanov/riemann_mod.f90',alias='riemann_{0}_mod'.format(d))

# mk.add(equations + '/two_point/ec/two_point_flux_mod.f90')

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
