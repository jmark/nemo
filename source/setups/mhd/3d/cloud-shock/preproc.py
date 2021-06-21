#!/usr/bin/env python

import sys,os

home = os.environ['NEMO_HOME']
setup = os.environ['NEMO_SETUP']
sys.path.append(os.path.join(home,'scripts'))

from nemo import Makefile

flags = dict(
    PP_N_NODES                  = 8,
    PP_KERNEL_DGFV_QMAX         = 2,
    PP_MESH_PERIODIC            = 0,
    DO_PROFILING                = 0,
)

mk = Makefile(home,flags)

#kernel = home + '/source/kernel/dgfv/3d'
kernel = home + '/source/kernel/mcfv/3d'

#timedisc = home + '/source/timedisc/rungekutta/ssp/4-5'
timedisc = home + '/source/timedisc/rungekutta/ralston'

equations = home + '/source/equations/mhd'

mk.add(equations + '/riemann/rusanov/riemann_mod.f90',alias='riemann_inner_mod')

for d in 'north south west east front back'.split():
    mk.add(equations + '/riemann/roe/riemann_mod.f90',alias='riemann_{0}_mod'.format(d))

mk.add(home + '/source/share_mod.f90')
mk.add(home + '/source/driver_mod.f90')
mk.add(home + '/source/runtime_mod.f90')

mk.add(home + '/source/utils/*.f90')
mk.add(home + '/source/mesh/**.f90')

mk.add(home + '/source/checkpoint/flash/*.f90')

mk.add(equations + '/*.f90')
mk.add(kernel + '/*.f90')
mk.add(timedisc + '/*.f90')

mk.add(setup + '/../_common_/*_mod.f90')
mk.add(setup + '/*_mod.f90')

mk.generate()
