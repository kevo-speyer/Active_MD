!
! Dipolar correction to the Ewald approach for coulombic interactions 
!
! Used when the system has slab geometry with walls in z=0 and z=D planes 
! Taken from Smit & Frenkel, chapter 12.4 and originally from 
! Yeh and Berkowitz:  

subroutine dipolar_correction()
#include 'control_simulation.h'
    use commons
    implicit none 
    real (kind=8) :: M_z
    logical, parameter :: debug=.false.

! Dipolar moment in Z

#if SYSTEM ==2 || SYSTEM ==3 


        M_z = 0.
        do i_part = 1,n_mon_tot
            M_z = M_z + q(i_part)*( r0(3,i_part)  - lz_half )            
        end do

! Dipolar force and energy

        do i_part = 1,n_mon_tot

            force(3,i_part) = force(3,i_part) - 4.*pi*q(i_part)*M_z*inv_V

        end do

        v_coul = v_coul + twopi*M_z*M_z*inv_V
        
        if(debug) then
            print '(a,3f16.5)' ,"[dipolar_correction] M_z,u_dipo,forze(q=+1) = ",M_z,  & 
                                        twopi*M_z*M_z*inv_V , -4.*pi*M_z*inv_V
        end if
#endif

end subroutine dipolar_correction
