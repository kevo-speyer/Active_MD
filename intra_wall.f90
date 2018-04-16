      subroutine intra_wall
      use commons ; implicit none
      real (kind=8) ::  dummy
      real (kind=8)          :: v_wall_intra
!
      v_wall_intra = 0.
!
!---  coupling to equilibrium sites
!obsol      if(f_wall_fix.eq.0) then
!obsol          do i_dim = 1,n_dim
!obsol              do i_part = n_mon_tot+1,n_mon_tot+n_wall/2
!obsol                  !----       top wall
!obsol                  r_dummy = k_spring_wall*                                        &
!obsol                      &      (r0(i_dim,i_part)-r_wall_equi(i_dim,i_part-n_mon_tot))
!obsol                  force(i_dim,i_part) = force(i_dim,i_part)-r_dummy
!obsol                  force_twall(i_dim) = force_twall(i_dim)+r_dummy
!obsol                  v_wall_intra = v_wall_intra + r_dummy**2
!obsol                  !----       bottom wall
!obsol                  j_part = i_part + n_wall/2
!obsol                  r_dummy = k_spring_wall*                                        &
!obsol                      &      (r0(i_dim,j_part)-r_wall_equi(i_dim,j_part-n_mon_tot))
!obsol                  force(i_dim,j_part) = force(i_dim,j_part)-r_dummy
!obsol                  v_wall_intra = v_wall_intra + r_dummy**2
!obsol              end do
!obsol          end do
!       else  !claudio: this is the case, for the usual use of the program

           do i_dim = 1,n_dim
               do i_part = n_mon_tot+1-n_wall,n_mon_tot-n_wall/2
                   force_twall(i_dim) = force_twall(i_dim)  + force(i_dim,i_part)
               end do
               ftw(i_dim)=force_twall(i_dim)
           end do

!      end if
!
      v_wall_intra = v_wall_intra/(2*k_spring_wall)
      v_wall_wall = v_wall_wall + v_wall_intra
!
#ifdef SHEARED
     do i_dim = 1,n_dim
       if(f_twall(i_dim).eq.1) then
!-----     constant or zero velocity of equilibrium positions
             f0_spring_1(i_dim) = f0_spring_1(i_dim) + force_twall(i_dim)
             force_twall(i_dim) = 0.
        end if
       end do
#endif
! WARN: compat with original program broken

! OBSOLETE !---- coupling to external spring
! OBSOLETE       do i_dim = 1,n_dim
! OBSOLETE        if(f_twall(i_dim).eq.1) then
! OBSOLETE !-----     constant or zero velocity of equilibrium positions
! OBSOLETE              f0_spring_1(i_dim) = f0_spring_1(i_dim) + force_twall(i_dim)
! OBSOLETE              force_twall(i_dim) = 0.
! OBSOLETE        else if(f_twall(i_dim).eq.2) then
! OBSOLETE !-----     add external force
! OBSOLETE            dummy = va_spring_twall(i_dim)
! OBSOLETE            va_spring_twall(i_dim) = va_spring_twall(i_dim) +               &
! OBSOLETE      &     dt * s_force_grad(i_dim) * k_spring_twall(i_dim)
! OBSOLETE            if(dummy*va_spring_twall(i_dim).lt.(0.d0))                      &
! OBSOLETE      &        s_force_grad(i_dim) = s_force_grad(i_dim)/4
! OBSOLETE               force_twall(i_dim) = force_twall(i_dim) +                       &
! OBSOLETE      &        n_wall/2 * va_spring_twall(i_dim)
! OBSOLETE            else if(f_twall(i_dim).eq.3) then
! OBSOLETE !-----     add force from external spring
! OBSOLETE              force_twall(i_dim) = force_twall(i_dim) + n_wall/2 *            &
! OBSOLETE      &       k_spring_twall(i_dim)*(r0_spring_twall(i_dim)-r0_twall(i_dim))
! OBSOLETE              f0_spring_1(i_dim) = f0_spring_1(i_dim) + n_wall/2 *            &
! OBSOLETE      &       k_spring_twall(i_dim)*(r0_spring_twall(i_dim)-r0_twall(i_dim))
! OBSOLETE         end if
! OBSOLETE        end do
!
      end subroutine intra_wall
