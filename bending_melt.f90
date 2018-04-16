! Routine that allows to add rigidity to the melt chains. CAlculates bending
! Forces and Energy.
! WARNING Bending constant(k_bend) is the same for grafted and melt chains, and
! is read from mfa_input file.
subroutine bending_melt(mode_bend_melt)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_bend_melt
integer ::i,j,l,beg
real(kind=8) ::r_prev(3),r_next(3),dot_pr=0,r_prev_2=0,r_next_2=0, cos_alpha,dir_next(3)=0,dir_prev(3)=0,dir_prev_2=0, dir_next_2=0, F_bend(6),F_mod,delta_alpha
#ifdef BENDING_MELT
select case (mode_bend_melt)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with bending stiffnes in melt"
print*," * Bending elastic constant for melt ",k_bend
print*," * Equilubrium angle for melt chains: ",alpha_eq
print*,""

case(1)
beg=n_mon*n_chain
v_bend_melt=0  !Reset Bending Potential Energy
Do l=0,n_chain_d-1 !n_chain , loop over melt chains
    r_next=r0(:,beg+l*n_mon_d+2)-r0(:,beg+l*n_mon_d+1)
    !This lines below correct for perioduc boundry conditions
    r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
    r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
   ! #if SYMMETRY == 1
   ! r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))
   ! #endif
    !End correction fo periodic boundry conditions
    Do i=2,n_mon_d-1 !n_mon, loop over particles in chain
        !Reset dummy variables
        dot_pr=0
        r_prev_2=0
        r_next_2=0
        dir_prev_2=0
        dir_next_2=0
        !Calculation of bending forces for the neighbors of i
        r_prev=r_next
        r_next=r0(:,beg+l*n_mon_d+i+1)-r0(:,l*n_mon_d+i)
        !This lines below correct for perioduc boundry conditions
        r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
        r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
        !#if SYMMETRY == 1
        !r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))
        !#endif
        !End correction fo periodic boundry conditions
        Do j=1,3
            dot_pr=dot_pr+r_next(j)*r_prev(j)
            r_prev_2=r_prev_2+r_prev(j)*r_prev(j)
            r_next_2=r_next_2+r_next(j)*r_next(j)
        End do
        cos_alpha=dot_pr/sqrt(r_prev_2*r_next_2)
        dir_next=r_prev*r_next_2-r_next*dot_pr
        dir_prev=-r_next*r_prev_2+r_prev*dot_pr
        Do j=1,3
            dir_prev_2=dir_prev_2+dir_prev(j)*dir_prev(j)
            dir_next_2=dir_next_2+dir_next(j)*dir_next(j)
        End do
        if(dir_next_2.lt.0.0001) then !if alpha is small make force 0
            dir_next=dir_next*0
        else
        !Below I divide by |r_next|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_next|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_next=dir_next/sqrt(dir_next_2*r_next_2)
        end if
        if(dir_prev_2.lt.0.0001) then !if alpha is small make force 0
            dir_prev=dir_prev*0
        else
        !Below I divide by |r_prev|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_prev|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_prev=dir_prev/sqrt(dir_prev_2*r_prev_2)
        end if
        delta_alpha=acos(cos_alpha)-alpha_eq
        F_mod=k_bend*delta_alpha
        v_bend_melt=v_bend_melt+.5*F_mod*delta_alpha
        Do j=1,3
            F_bend(j)=F_mod*dir_prev(j)
            F_bend(3+j)=F_mod*dir_next(j)
        End do
        Do j=1,3
            force(j,beg+l*n_mon_d+i-1) = force(j,beg+l*n_mon_d+i-1) + F_bend(j)
            force(j,beg+l*n_mon_d+i)   = force(j,beg+l*n_mon_d+i) - F_bend(j) - F_bend(3+j)
            force(j,beg+l*n_mon_d+i+1) = force(j,beg+l*n_mon_d+i+1) + F_bend(3+j)
        End do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  DEBUGGING TO TEST TO CHECK IF THE ROUTINE IS CALCULATING FORCES AND ENERGY
        !  CORRECTLY. 
        !print*, r0, "positions"
        !print*, k_bend, "bending constant"
        !print*, F_bend, "bending FORCE"
        !print*, v_bend_melt, "bending ENERGY"
        ! ERASE AFTER CHECKING
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
    End Do
End do
end select

#endif /*close BENDING_MELT*/
end subroutine bending_melt
