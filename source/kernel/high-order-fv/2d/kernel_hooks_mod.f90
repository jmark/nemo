module kernel_hooks_mod

use constants_mod

contains

subroutine hook_pre_rk_timestep_300

    use share_mod, only: mesh

    use mesh_utils_mod, only: mesh_iterate_quads
    use mesh_utils_mod, only: mesh_iterate_sides
    use mesh_const_mod, only: MESH_ITERATE_WITH_GHOSTS

    use kernel_flux_mod, only: kernel_fill_halo

    ! call mesh_iterate_quads(mesh,kernel_fill_facevars,flags=MESH_ITERATE_WITH_GHOSTS)
    call mesh_iterate_quads(mesh,kernel_fill_halo)

end subroutine

end module
