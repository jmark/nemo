module timedisc_evolve_mod

use constants_mod

# if DO_PROFILING
integer, save :: prof_rkstep
# endif

contains

subroutine hook_init_399

# if DO_PROFILING
    use profiler_mod, only: profiler_init
    prof_rkstep = profiler_init('rkstep')
# endif

end subroutine

subroutine hook_evolve_500

    use share_mod, only: mesh
    use share_mod, only: rt => runtime
    use mesh_utils_mod, only: mesh_iterate_quads

    use hook_pre_timestep_mod, only: hook_pre_timestep

    use hook_pre_rk_timestep_mod, only: hook_pre_rk_timestep
    use hook_mid_rk_timestep_mod, only: hook_mid_rk_timestep
    use hook_post_rk_timestep_mod, only: hook_post_rk_timestep

    use hook_post_timestep_mod, only: hook_post_timestep

    use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_INTEGER
    use status_mod, only: raise_error

    use timedisc_const_mod, only: nstages
    use timedisc_const_mod, only: dt_coeffs

# if DO_PROFILING
    use profiler_mod, only: profiler_start, profiler_stop
# endif

    logical :: ok
    integer :: stage

    integer :: locstat,glostat !! local and global status code

    call hook_pre_timestep()

    do stage = 1,nstages
        rt%rkstage    = stage
        rt%rktimestep = rt%timestep * dt_coeffs(stage)
        rt%rksimtime  = rt%simtime  + rt%rktimestep

        call hook_pre_rk_timestep()
# if DO_PROFILING
        call profiler_start(prof_rkstep)
# endif
        call hook_mid_rk_timestep
# if DO_PROFILING
        call profiler_stop(prof_rkstep)
# endif
        call hook_post_rk_timestep

        ! call raise_error('Debugging.')

        if (rt%mpi%size > 1) then
            locstat = merge(0,1,rt%stat%ok)
            glostat = 0
            call MPI_Allreduce(locstat,glostat,1,MPI_INTEGER,MPI_SUM,rt%mpi%comm)

            if (glostat > 0) then
                if (rt%stat%ok) call raise_error('An error on another MPI rank occurred.')
                return
            end if
        else
            if (.not.rt%stat%ok) return
        end if
    end do

    call hook_post_timestep

end subroutine

end module
