module share_mod

    use mesh_types_mod, only: mesh_t
    use runtime_mod, only: runtime_t
    use traceback_mod, only: traceback_t
# if PP_CHECKPOINT
    use checkpoint_types_mod, only: checkpoint_t
# endif

    type(runtime_t), save, target   :: runtime
    type(traceback_t), save, target :: traceback

    !! Pointer to current active mesh.
    type(mesh_t), save, pointer     :: mesh

    !! Necessary for AMR. One is the old mesh
    !! while the other is the new mesh.
    type(mesh_t), save, target :: meshes(2)
    integer, save :: mesh_active = 1

# if PP_CHECKPOINT
    type(checkpoint_t), save, target :: chkpt
# endif

end module
