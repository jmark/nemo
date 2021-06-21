module FMM3D_mod

use constants_mod
     
interface
    subroutine lfmm3d_t_c_p(eps,nsources,sources,charges,ntargets,targets,potential,plummer)

        import dp

        real(dp), intent(in)    :: eps
        integer, intent(in)     :: nsources,ntargets
      
        real(dp), intent(in)    :: sources(3,nsources)
        real(dp), intent(in)    :: charges(  nsources)

        real(dp), intent(in)    :: targets(3,ntargets)
        real(dp), intent(out)   :: potential(ntargets)

        real(dp), intent(in)    :: plummer

    end subroutine
end interface

interface
    subroutine lfmm3d_t_c_g(eps,nsources,sources,charges,ntargets,targets,potential,gradients,plummer)

        import dp

        real(dp), intent(in)    :: eps
        integer, intent(in)     :: nsources,ntargets
      
        real(dp), intent(in)    :: sources(3,nsources)
        real(dp), intent(in)    :: charges(  nsources)

        real(dp), intent(in)    :: targets(3,  ntargets)
        real(dp), intent(out)   :: potential(  ntargets)
        real(dp), intent(out)   :: gradients(3,ntargets)

        real(dp), intent(in)    :: plummer

    end subroutine
end interface

interface
    subroutine lfmm3d_st_c_g(eps,nsources,sources,charges,potl,grad,ntargets,targets,potltarg,gradtarg,plummer)

        import dp

        real(dp), intent(in)    :: eps
        integer, intent(in)     :: nsources,ntargets
      
        real(dp), intent(in)    :: sources(3,nsources)
        real(dp), intent(in)    :: charges(  nsources)

        real(dp), intent(out)   :: potl(  ntargets)
        real(dp), intent(out)   :: grad(3,ntargets)

        real(dp), intent(in)    :: targets(3, ntargets)
        real(dp), intent(out)   :: potltarg(  ntargets)
        real(dp), intent(out)   :: gradtarg(3,ntargets)

        real(dp), intent(in)    :: plummer

    end subroutine
end interface

end module
