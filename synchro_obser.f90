subroutine synchro_obser(mode)
#include "control_simulation.h"
    use commons
    implicit none
    INTEGER, intent(in) :: mode
    INTEGER :: i_ch

    select case(mode)

    case(1)

        write(56,'(I12.4)', advance='no') i_time
#ifdef ACTIVE_BRUSH
        write(56,'(ES12.4)', advance='no') cos_mem(n_chain)
        write(56,'(ES12.4)', advance='no') k0(1+(n_chain-1)*(n_mon-1))
#else
        r_rel = r0(:, 2+(n_chain-1)*n_mon)- r0(:, 1+(n_chain-1)*n_mon)  
        cos_th = r_rel(1)/SQRT(DOT_PRODUCT(r_rel,r_rel))
        write(56,'(ES12.4)', advance='no') cos_th
        write(56,'(ES12.4)', advance='no') k_bend  ! Writing useless column for script compatibility
#endif
        ! do i=1,3
        !     write(56,'(ES12.4)', advance='no') f_old(i)-f_new(i)        
        ! end do
        write(56,'(ES12.4)', advance='no') pot
        write(56,*)

    case(2)

        write(56,'(I12.4)', advance='no') i_time
        do i_ch=1,n_chain
            r_rel = r0(:, n_mon*i_ch)- r0(:, 1+(n_chain-1)*n_mon)  
            cos_th = r_rel(1)/SQRT(DOT_PRODUCT(r_rel,r_rel))
            write(56,'(ES12.4)', advance='no') cos_th
        enddo
        write(56,*)
    end select
end subroutine synchro_obser
