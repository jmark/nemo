module mesh_hooks_mod

use constants_mod

contains

subroutine hook_init_500

    use share_mod, only: mesh, meshes, mesh_active
    use mesh_utils_mod, only: mesh_init

    mesh => meshes(1)
    mesh_active = 1

    call mesh_init(mesh)

end subroutine

subroutine hook_nuke_800

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_nuke

    call mesh_nuke(mesh)

end subroutine

subroutine hook_pre_rk_timestep_200

    use share_mod, only: mesh
    use mesh_utils_mod, only: mesh_exchange

    call mesh_exchange(mesh)

end subroutine

end module
