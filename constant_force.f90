subroutine constant_force
#include 'control_simulation.h' 
! This routine adds an external constant force to the  monomers.
! Gravity or force in x. If the masses are different, the force should be
! thought as m*g force, where g is read from the force parameter in
! system_input
! The force magnitud is constant and read form system_input
! * the force is added to the particles with a_type = 3 
    use commons
!use ziggurat, only: rnor,uni
!use util ! debug
    implicit none
    logical, parameter :: debug=.false.
    
    
    do i_part = 1, n_mon_tot
          force(:,i_part) = force(:,i_part) + ext_force(:,i_part)
    end do

end subroutine constant_force
