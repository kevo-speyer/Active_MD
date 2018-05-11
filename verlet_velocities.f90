subroutine verlet_velocities()
! * Updates velocities, forces and accels using Velocity verlet  
! * Updates wall positions
! in the DPD scheme
use commons
#include "control_simulation.h" 
  implicit none
  real (kind=8) ::  T_inst, T_dummy
  real (kind=8) ::  T_inst_mean = 0.
  integer :: i_head
  logical, parameter :: debug = .false.


#if SYMMETRY == 0
#       if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3 /* all the systems with fixed brush heads */
#        ifndef FREE_HEADS
!   ---  Zero all for the heads
                f_on_heads(:,:) =  force(:,1:part_init_d:n_mon)
#         endif
#       endif
#endif

! Update accelerations with the new force values 

        do i_part = 1,n_mon_tot
            a(1,i_part) = force(1,i_part)*inv_mass(i_part)
#ifdef BIDIMENSIONAL
            a(2,i_part) = 0.0d0 
#else 
            a(2,i_part) = force(2,i_part)*inv_mass(i_part)
#endif
            a(3,i_part) = force(3,i_part)*inv_mass(i_part)
!        end do

! Update velocities 

!        do i_part = 1 , n_mon_tot
!            v(1,i_part) = v_half(1,i_part) + 0.5*dt*a(1,i_part)       ! old force(1,i_part)*inv_mass(i_part)
!            v(2,i_part) = v_half(2,i_part) + 0.5*dt*a(2,i_part)       ! old force(2,i_part)*inv_mass(i_part)
!            v(3,i_part) = v_half(3,i_part) + 0.5*dt*a(3,i_part)       ! old force(3,i_part)*inv_mass(i_part)
            v(1,i_part) = v(1,i_part) + 0.5*dt*a(1,i_part)    
#ifdef BIDIMENSIONAL
            v(2,i_part) = 0.0d0
#else 
            v(2,i_part) = v(2,i_part) + 0.5*dt*a(2,i_part)    
            !print *, v(2,i_part), a(2,i_part)
#endif
            v(3,i_part) = v(3,i_part) + 0.5*dt*a(3,i_part)    
        end do


! NOTE: if DPD_VV is defined the energies and temp are calculated afterwards in new_dpd_fd.f90
#ifdef SHEARED
    if(f_twall(i_dim).eq.9) then
        do i_dim=1,n_dim
            appl_vel_tw(i_dim)=appl_vel_tw(i_dim)+0.5*dt_2*ftw(i_dim)/mass(n_part)
            appl_vel_bw(i_dim)=appl_vel_bw(i_dim)+0.5*dt_2*fbw(i_dim)/mass(n_part)
        end do
    end if
#endif
#ifndef DPD_VV


#ifdef SHEARED
    do i_dim=1,n_dim
        if(f_twall(i_dim).eq.9) then
            appl_acc_tw(i_dim) = 0.5*dt_2*ftw(i_dim)/mass(n_part)
            appl_acc_bw(i_dim) = 0.5*dt_2*fbw(i_dim)/mass(n_part)
        end if
    end do
#endif

#ifdef RELAX
! If doing RELAXATION WITHOUT DYNAMICS 

! Velocities to 0 in each time step  and relax in the direction of  the forces !!!

           v(:,:) =0.

! RELAX Constraint the maximum force to 10.

            do i_part = 1,n_mon_tot

                if ( force(1,i_part) > 10. )   force(1,i_part) = 10.
                if ( force(2,i_part) > 10. )   force(2,i_part) = 10.
                if ( force(3,i_part) > 10. )   force(3,i_part) = 10.

                if ( force(1,i_part) < -10. )  force(1,i_part) = -10.
                if ( force(2,i_part) < -10. )  force(2,i_part) = -10.
                if ( force(3,i_part) < -10. )  force(3,i_part) = -10.

            end do
#endif /* ifdef RELAX */



!!KEVIN 10/2013 Droplet compatibility added (WALL ==3 and SYSTEM == 1) 
!!Droplet not compatible with moving wall

#if SYMMETRY == 0 /* channel geometry */
#   ifndef POISEUILLE /* if not defined poiseuille */
#      if WALL  == 1 || WALL == 2 || WALL == 3 /* if using or explicit implicit wall with fixed brush heads */
#           if SYSTEM == 0 || SYSTEM == 3 || SYSTEM == 1
!
!       ----- Wall velocities,accels and forces for brush's head 
!
        do i_chain = 1 ,n_chain/2
                i_head = 1+n_mon*(i_chain-1)
                v(1:2,i_head) = va_spring_twall(1:2)
                v(3,i_head) = 0.
                force(:,i_head) = 0.! zero force
                a(:,i_head) = 0. ! zero accels
        end do
        do i_chain = n_chain/2+1 ,n_chain
                i_head = 1+n_mon*(i_chain-1)
                v(1:2,i_head) = -va_spring_twall(1:2)
                v(3,i_head) = 0.
                force(:,i_head) = 0.! zero force
                a(:,i_head) = 0. ! zero accels
        end do
#           endif
#       endif
#   else /* if Poiseuille is defined */
!      Fix the head beads at the original position 
        do i_chain = 1 ,n_chain
                i_head = 1+n_mon*(i_chain-1)
                v(:,i_head) = 0. ! 0 the velocity of the heads
                force(:,i_head) = 0.! zero force
                a(:,i_head) = 0. ! zero accels
        end do
#   endif /* indef poiseuille*/
# endif /* symmetry */

#ifndef RESPA

!NOW THIS IS DONE IN verlet_velocities_long.f90 

!   ---- Measurement of some quantities 


!---  Kinetic energies and temperatures 


! Instantaneous temperature

      T_inst = 0.
      do i_part = 1,n_mon_tot
          T_dummy =   mass(i_part)*( v(1,i_part)**2+ v(2,i_part)**2 + v(3,i_part)**2 )
          T_inst = T_inst + T_dummy
      end do
      t_fluid = T_inst
  
#   if     WALL == 1     /* explicit wall */
        t_wall = 0.
        do i_part = n_mon_tot+1,n_part
            do i_dim= 1,n_dim
                t_wall = t_wall + mass(i_part)*v(i_dim,i_part)**2
            end do
        end do
      t_wall = 0.5*t_wall ! *inv_dt_2
#   endif        


! *** Write out Instantanous T and energies
!

      t_fluid = 0.5*t_fluid ! *inv_dt_2
!
      v_total = v_wall_wall + v_fluid_wall + v_fluid_fluid + v_intra_molec

#ifdef BENDING
        v_total = v_total + v_bend !Add bending energy  
#endif

#ifdef ORIENTATION
        v_total = v_total + v_or !Add bending orientation energy
#endif 

#ifdef BENDING_MELT
      v_total = v_total + v_bend_melt !Add melt bending energy
#endif

#ifdef SPRING_ARRAY
      v_total = v_total + v_array  ! Add spring array energy
#endif

!        print *,"v_fluid_fluid=",v_fluid_fluid/dble(n_part) !; stop

!
! Update wall velocity 

      if(f_twall(n_dim).eq.2) then
          v_total=v_total-r0_twall(n_dim)*va_spring_twall(n_dim)*0.5
      end if
!
      t_total = t_wall + t_fluid
      e_total = v_total + t_total
      
!      
!!  Temperature and energy calculation and writing out 
!

   time_ave_count_2 = time_ave_count_2 + 1

   if(time_ave_count_2.eq.10) then
       time_ave_count_2 = 0

#ifdef SHEARED
!Claudio:  I assume that this out of sheared is wrong

        r0_twall_1(:) = r0_twall_1(:) / dble(n_time_ave)
        f0_spring_1(:) = f0_spring_1(:)/ dble(n_time_ave)
        r0_twall_1(:) = 0.
        f0_spring_1(:) = 0.

#endif        

#if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 4
        write(61,'(i7,3g17.5)') i_time,v_fluid_fluid,v_intra_molec,v_fluid_wall
        write(60,'(i11,3g17.5)') i_time,e_total*inv_N,v_total*inv_N, t_total*inv_N
#elif SYSTEM == 2  || SYSTEM == 3
        v_total = v_total + v_coul
        write(60,'(i11,4g17.5)') i_time,e_total*inv_N,v_total*inv_N, t_total*inv_N,v_coul*inv_N
#endif

    end if
#endif   /* not DPD_VV */ 


!
!   *** DEBUGGING STUFF *** 
!
if(debug) then

print '(a,4e15.5)',"[Venergies]" , v_wall_wall , v_fluid_wall , v_fluid_fluid , v_intra_molec
print*,"before_vel^2= ",i_time,sum(v(:,1:n_mon_tot)**2, dim=1)*inv_N
print*,"V_tot_Ekin",i_time,v_total*inv_N,t_total*inv_N
print*,"heads_vel= ",sum(v(:,1:n_mon*n_chain:n_mon), dim=1)
print*,"accel^2= ",i_time,sum(a(:,1:n_mon_tot)**2, dim=1)/real(n_mon_tot)
        write(77,'(i7,3f17.5)') i_time,v_total*inv_N, e_total*inv_N,t_total*inv_N
end if

#endif /* ifndef RESPA*/
end subroutine verlet_velocities
