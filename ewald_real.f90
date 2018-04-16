    subroutine ewald_real()
#include 'control_simulation.h' 
! This routine computes, ewald sum implementation of the Coulomb interactions
! in the real space
!
   use commons
   implicit none 
!   real(kind=8),intent(in) :: r_2
!   integer, intent(in) :: i_part,j_part
!   real(kind=8),intent(in) :: dr(3)
    real(kind=8)            :: r,inv_r, f_cou_r(3),f_scal,qq
    logical, parameter :: debug=.true.

#if SYSTEM == 2 || SYSTEM == 3
    r = sqrt(r_2)
    inv_r = 1./r 

    qq = q(i_part)*q(j_part)
      
    v_coul = v_coul + qq*erfc(sqrt_alpha*r)*inv_r
     
    f_scal = qq*sqrt_alpha*inv_sqrt_pi*2.*exp(-alpha_c*r_2)*inv_r + qq*erfc(sqrt_alpha*r)*inv_r*inv_r

!NOTE: the factor 1/2 is omitted because we visit each pair of atoms only once. 

    f_cou_r(:) = f_scal*delta_r(:)*inv_r 

! Add to the total force 

       force(:,i_part) = force(:,i_part)  + f_cou_r(:)
#     if BIN_TYPE == 0       
       force(:,j_part) = force(:,j_part)  - f_cou_r(:)
#     endif
    
    if(debug) then
!    print *,"[ewald_real] ",i_time, v_coul 
    end if

#endif

    end subroutine ewald_real
