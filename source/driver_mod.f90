module driver_mod

contains

subroutine mainloop()

    use constants_mod, only: dp,TOL,N_NODES,N_DIMS
    use share_mod, only: mesh, rt => runtime 

# if PP_CHECKPOINT
    use hook_checkpoint_mod, only: hook_checkpoint
# endif
    use hook_calc_timestep_mod, only: hook_calc_timestep
    use hook_analyze_mod, only: hook_analyze
    use hook_evolve_mod, only: hook_evolve

    use timer_mod, only: timer_t
    use timer_mod, only: timer_start
    use timer_mod, only: timer_stop
    use timer_mod, only: timer_read

    use setup_config_mod, only: init_time
    use setup_config_mod, only: stop_time
    use setup_config_mod, only: dt_analyze
    use setup_config_mod, only: dt_checkpoint
    use setup_config_mod, only: chkpt_init_time

    use timedisc_const_mod, only: nstages

# if DO_PROFILING
    use setup_config_mod, only: prof_report_rate
    use profiler_mod, only: profiler_init,profiler_start, profiler_stop, profiler_report

    integer :: prof_cycle,prof_analyze,prof_checkpoint,prof_timestep,prof_beat,prof_pid
# endif

    integer :: heartbeat    = 0
    integer :: lastdump     = 0
    integer :: nsteps       = 1

    type(timer_t)   :: timer_pid !! performance index
    real(dp)        :: perfidx

# if DO_PROFILING
    prof_beat       = 0
    prof_cycle      = profiler_init('cycle')
    prof_analyze    = profiler_init('analyze')
    prof_timestep   = profiler_init('timestep')
    prof_checkpoint = profiler_init('checkpoint')
# endif

    rt%simtime      = init_time
    rt%timestep     = 0.0_dp
    rt%timesteps    = 0
    rt%nextdump     = rt%simtime

    if (rt%mpi%root) write(*,1000) 'info    ', '#io', '#dt', 'simtime ', 'timestep', 'perf. index'
    if (rt%mpi%root) write(*,*)

    call timer_start(timer_pid)

    do while (rt%stat%ok)

# if DO_PROFILING
        call profiler_start(prof_cycle)
        call profiler_start(prof_timestep)
# endif
        call hook_calc_timestep
# if DO_PROFILING
        call profiler_stop(prof_timestep)
# endif

# if DO_PROFILING
        call profiler_start(prof_analyze)
# endif
        call hook_analyze
# if DO_PROFILING
        call profiler_stop(prof_analyze)
# endif

        if (abs(rt%simtime - rt%nextdump) < 10*TOL) then
            call timer_stop(timer_pid)
            perfidx = timer_read(timer_pid) * REAL(rt%ncpus,dp) / REAL(mesh%nglobal * N_NODES**N_DIMS * heartbeat * nstages,dp) * 1e6_dp
            if (rt%mpi%root) write(*,1001) rt%chkptcount, rt%timesteps, rt%simtime, rt%maxtimestep, perfidx

# if DO_PROFILING
            call profiler_start(prof_checkpoint)
# endif

# if PP_CHECKPOINT
            call hook_checkpoint()
# endif
# if DO_PROFILING
            call profiler_stop(prof_checkpoint)
# endif

            heartbeat    = 0
            lastdump     = 0

            rt%nextdump = min(stop_time,rt%simtime + dt_checkpoint)

            call timer_start(timer_pid)
        end if

        if (abs(rt%simtime - stop_time) < 100*TOL) exit

        rt%timestep = min(&
            rt%maxtimestep, &
            dt_checkpoint, &
            dt_analyze, &
            abs(stop_time - rt%maxtimestep), &
            abs(rt%nextdump - rt%simtime))

        call hook_evolve()

# if DO_PROFILING
        call profiler_stop(prof_cycle)

        prof_beat = prof_beat + 1
        if (prof_beat >= prof_report_rate) then
            call profiler_report
            prof_beat = 0
        end if
# endif

        if (.not.rt%stat%ok) then
# if PP_CHECKPOINT
            if (rt%mpi%root) write(*,1001) rt%chkptcount, rt%timesteps, rt%simtime, rt%maxtimestep
            call hook_checkpoint()
# endif
            exit
        end if

        rt%simtime      = rt%simtime    + rt%timestep
        rt%timesteps    = rt%timesteps  + 1
        heartbeat       = heartbeat     + 1
        lastdump        = lastdump      + 1
        nsteps          = nsteps        + 1

        if (heartbeat >= 64) then
            call timer_stop(timer_pid)
            perfidx = timer_read(timer_pid) * REAL(rt%ncpus,dp) / REAL(mesh%nglobal * N_NODES**N_DIMS * heartbeat * nstages,dp) * 1e6_dp

            if (rt%mpi%root) then
                write(*,1002) rt%chkptcount-1, rt%timesteps, rt%simtime, rt%maxtimestep, perfidx
            end if
            heartbeat = 0
            call timer_start(timer_pid)
        end if
    end do

1000 format (a12, a7, 4x, a8, 2x, a11, 2x, a11, 4x, a11)
1001 format ('[CHECKPOINT]', i6, 4x, i8, 2x, ES12.3, 2x, ES12.4, 4x,F8.4)
1002 format ('[ HEARTBEAT]', i6, 4x, i8, 2x, ES12.3, 2x, ES12.4, 4x,F8.4)

end subroutine

end module
