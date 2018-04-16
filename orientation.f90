! Routine that allows to add rigidity to the 2 bead of the grafted chains, in
! order to induce an orientation angle
subroutine orientation(mode_or)
#include "control_simulation.h"
use commons
implicit none
integer ,intent(in) :: mode_or 
integer ::i,j,l
real(kind=8) ::r_prev(3),r_next(3),dot_pr=0,r_prev_2=0,r_next_2=0,cos_alpha,dir_next(3)=0,dir_prev(3)=0,dir_prev_2=0, dir_next_2=0,F_mod,delta_alpha,r_ghost(3)
#ifdef ORIENTATION
select case (mode_or)

case(0)  ! Init  variables 
#ifdef ACTIVE_BRUSH
    ALLOCATE(k0( (n_mon-1)*n_chain) )
    k0 = k_bend 
#else
    k_or = k_bend!sets orientation stiffnes
#endif
    alpha_or=0
    print*,""
    print*," * Simulation with orientation stiffnes"
    print*," * Orientation elastic constant ",k_or
    print*," * Orientation angle:",alpha_or
    print*," * Force on heads correctly calculated"
    print*," * Valid fort.80"
     

case(1)  ! orientation rigidity calculation
    v_or = 0  !Reset orientation Potential Energy

    Do l=0,n_chain-1 !n_chain , loop over chains
        !sets r_ghost for orientation stiffnes calculation
        dot_pr=0
        r_prev_2=0
        r_next_2=0
        dir_prev_2=0
        dir_next_2=0 

        r_ghost=r0(:,l*n_mon+1)
        if(r0(3,l*n_mon+1).gt.boundary(3)/2) then
            r_ghost(3) = r_ghost(3) + 1.0
        else
            r_ghost(3) = r_ghost(3) - 1.0     
        end if
        r_prev = r0(:,l*n_mon+1) - r_ghost
        r_next = r0(:,l*n_mon+2) - r0(:,l*n_mon+1)
        !This lines below correct for periodic boundry conditions
        r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
        r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
        r_prev(1) = r_prev(1) - boundary(1) * int(2.*r_prev(1)*inv_boundary(1))
        r_prev(2) = r_prev(2) - boundary(2) * int(2.*r_prev(2)*inv_boundary(2))

        Do j=1,3
            dot_pr=dot_pr+r_next(j)*r_prev(j)
            r_prev_2=r_prev_2+r_prev(j)*r_prev(j)
            r_next_2=r_next_2+r_next(j)*r_next(j)
        End do
        cos_alpha=dot_pr/sqrt(r_prev_2*r_next_2)
        dir_next=r_prev*r_next_2-r_next*dot_pr
        !dir_prev=-r_next*r_prev_2+r_prev*dot_pr
        Do j=1,3
            !dir_prev_2=dir_prev_2+dir_prev(j)*dir_prev(j)
            dir_next_2=dir_next_2+dir_next(j)*dir_next(j)
        End do
        if(dir_next_2.lt.0.0001) then !if alpha is small make force 0
            dir_next=dir_next*0
        else
        !Below I divide by |r_next|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_next|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_next=dir_next/sqrt(dir_next_2*r_next_2)
        end if
        !if(dir_prev_2.lt.0.0001) then !if alpha is small make force 0
        !    dir_prev=dir_prev*0
        !else
        !    dir_prev=dir_prev/sqrt(dir_prev_2)
        !end if
        delta_alpha=acos(cos_alpha)-alpha_or
        if(cos_alpha.ge.1.0d0) then ! changed 1.0 to 1.0d0
                delta_alpha=0.0-alpha_eq
        end if
#ifdef ACTIVE_BRUSH
        F_mod = k0(1+l*(n_mon-1))*delta_alpha
#else
        F_mod = k_or*delta_alpha
#endif
        v_or=v_or+.5*F_mod*delta_alpha
        Do j=1,3
            force(j,l*n_mon+1) = force(j,l*n_mon+1) - F_mod*dir_next(j)
            force(j,l*n_mon+2) = force(j,l*n_mon+2) + F_mod*dir_next(j)
        End do

    End do

end select

#endif /*close inlination*/
end subroutine orientation



