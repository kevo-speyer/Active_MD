!
! This routine maintains the geometry of the droplet during thermalization
! It does only a geometric checking and put the particles back inside de box  
!
subroutine geom_constraint() 
use commons
implicit none

#if SYSTEM == 1 /*droplets */

! WARN: the vapour chains must be excluded of the geometric reshaping 
!       of the droplet
! ver2: we do it reversing velocities
        do i_part = part_init_d+1 , n_mon_tot ! - n_mon_vap
            if (r0(1,i_part) < x_box_min .and. v(1,i_part) < 0. ) then
                v(1,i_part) = -v(1,i_part)
            end if
            if (r0(1,i_part) > x_box_max .and. v(1,i_part) > 0. ) then
                v(1,i_part) = -v(1,i_part)
            end if
            if (r0(3,i_part) > boundary_d(3) .and. v(3,i_part) > 0.) then
                v(3,i_part) = -v(3,i_part)
            end if
        end do
#endif /*SYSTEM == 1 droplets */

end subroutine geom_constraint

