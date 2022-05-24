module mesh_const_mod

use constants_mod

integer, parameter :: MESH_DIMS     = N_DIMS
integer, parameter :: MESH_NODES    = N_NODES
integer, parameter :: MESH_HALF     = 2**(MESH_DIMS-1)
integer, parameter :: MESH_FACES    = 2*MESH_DIMS
integer, parameter :: MESH_CHILDREN = 2**(MESH_DIMS)
integer, parameter :: MESH_CORNERS  = 2**(MESH_DIMS)
integer, parameter :: MESH_VERTS    = 2**(MESH_DIMS)

integer, parameter :: MESH_BOUNDARY_INNER = 0
integer, parameter :: MESH_BOUNDARY_OUTER = 1

integer(kind=1), parameter :: MESH_FLAG_NONE    = 0
integer(kind=1), parameter :: MESH_FLAG_REFINE  = 1
integer(kind=1), parameter :: MESH_FLAG_COARSEN = 2

integer, parameter :: MESH_ITERATE_WITH_GHOSTS          = 1
integer, parameter :: MESH_ITERATE_ONLY_INNER           = 2
integer, parameter :: MESH_ITERATE_ONLY_GHOSTS          = 4
integer, parameter :: MESH_ITERATE_ONLY_ZONE_BOUNDARY   = 8
integer, parameter :: MESH_ITERATE_ONLY_INNER_WITHOUT_MIRRORS = 16

#:def MAPPING(**mapping)
#:set VERTICES = dict(NWF=1,SWF=2,NEF=3,SEF=4,NWB=5,SWB=6,NEB=7,SEB=8)
$:','.join(mapping[k] for k in sorted(mapping,key=lambda v: VERTICES[v]))
#:enddef

integer, parameter :: map_faceid_to_direction(6) = (/X_DIR,X_DIR,Y_DIR,Y_DIR,Z_DIR,Z_DIR/)
integer, parameter :: map_faceid_to_opposite(6) = (/SOU,NOR,EAS,WES,BAC,FRO/)

# if PP_N_DIMS == 2

integer, parameter :: map_to_opposite_vertid(MESH_CORNERS) = &
    (/${MAPPING(NWF='SEF',SWF='NEF',NEF='SWF',SEF='NWF')}$/)

integer, parameter :: map_faceid_to_vertids(2,MESH_FACES) = reshape(&
    (/ NWF,NEF, SWF,SEF, NWF,SWF, NEF,SEF /),(/2,MESH_FACES/))

integer, parameter :: map_faceid_to_hanging_vertid_verts(2,MESH_FACES) = reshape(&
    (/ SWF,SEF, NWF,NEF, NEF,SEF, NWF,SWF /),(/2,MESH_FACES/))

integer, parameter :: map_faceid_to_hanging_vertid_quads(2,MESH_FACES) = reshape(&
    (/ NEF,NWF, SEF,SWF, SWF,NWF, SEF,NEF /),(/2,MESH_FACES/))

# elif PP_N_DIMS == 3

integer, parameter :: map_to_opposite_vertid(MESH_CORNERS) = &
    (/${MAPPING(NWF='SEB',SWF='NEB',NEF='SWB',SEF='NWB',NWB='SEF',SWB='NEF',NEB='SWF',SEB='NWF')}$/)

integer, parameter :: map_faceid_to_vertids(4,MESH_FACES) = reshape((/&
    NWF,NEF,NWB,NEB, SWF,SEF,SWB,SEB, NWF,SWF,NWB,SWB, NEF,SEF,NEB,SEB, NWF,SWF,NEF,SEF, NWB,SWB,NEB,SEB &
/),(/4,MESH_FACES/))

integer, parameter :: map_faceid_to_hanging_vertid_verts(4,MESH_FACES) = reshape((/&
    SWF,SEF,SWB,SEB, NWF,NEF,NWB,NEB, NEF,SEF,NEB,SEB, NWF,SWF,NWB,SWB, NWB,SWB,NEB,SEB, NWF,SWF,NEF,SEF &
/),(/4,MESH_FACES/))

integer, parameter :: map_faceid_to_hanging_vertid_quads(4,MESH_FACES) = reshape((/&
    NEB,NWB,NEF,NWF, SEB,SWB,SEF,SWF, SWB,NWB,SWF,NWF, SEB,NEB,SEF,NEF, SEF,NEF,SWF,NWF, SEB,NEB,SWB,NWB &
/),(/4,MESH_FACES/))

# else
# error: Unsupported dimension!
# endif

end module
