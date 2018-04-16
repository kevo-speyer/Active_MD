subroutine verlet_positions()
#include 'control_simulation.h'
!  This is the first step velocity verlet 
!  
! * Updates coordinates 
! * Velocity to t + 0.5*dt 
use commons
!use util, only: write_conf
!use util ! only for generating the histogram, remove when not debugging
!use ziggurat !gaussian random generator used
!use util
      implicit none
      real (kind=8), pointer :: r_head(:,:) !(from predict.f90)
      real (kind=8) :: rr_dummy(3)
      logical, parameter :: debug=.false.

!! NOTE: now the head velocities are correctly set up in the velocity upgrade, so that
!       the positions will be updated here in the usual VV procedure
! If not explicit wall atoms, the brushes heads are fixed in r0 init 
!
! Remember the head positions before updating  r0
! if(.not.f_explicit_wall.and.n_chain>0) then
!claudio 
! -- Remember initial coordinates for brushes heads

! r_head_old(:,:) = r0(:,1:n_mon*n_chain:n_mon)
! -- Take the heads coor which are not integrated
!     r_head => r0(:,1:n_chain*n_mon:n_mon)
! end if 
     
#       ifdef SHEARED

    if(i_time.eq.1) then
        do i_dim = 1, n_dim
            appl_vel(i_dim)=va_spring_twall(i_dim)

            if(f_twall(i_dim).eq.0.or.f_twall(i_dim).eq.9) then
                appl_vel(i_dim)=0.
                appl_vel_tw(i_dim)=0.
                appl_vel_bw(i_dim)=0.
            end if
    
        end do 
    end if
    
    do i_dim = 1, n_dim
        if(f_twall(i_dim).eq.0.and.i_time.eq.turn_time(i_dim)) then
            appl_vel(i_dim)=va_spring_twall(i_dim)
            write(*,*) "applying velocity at time step",i_time
        end if
        if(f_twall(i_dim).eq.2) then
            if(i_time.eq.turn_time(i_dim)) then
                appl_vel(i_dim)=- va_spring_twall(i_dim)
                write(*,*) "shear direction inverted at time step",i_time
            end if
        end if
        if(f_twall(i_dim).eq.3) then

            if(i_time.eq.turn_time(i_dim)/2) then
                appl_vel(i_dim)=- va_spring_twall(i_dim)
                write(*,*) "shear direction inverted at time step",i_time
            end if
            if(mod(i_time+turn_time(i_dim)/2,turn_time(i_dim)).eq.0.and.i_time.gt.turn_time(i_dim)/2) then
                appl_vel(i_dim)=-appl_vel(i_dim)
                write(*,*) "shear direction inverted at time step",i_time
            end if
        end if


        if(i_time.eq.turn_time(i_dim)) then
            if(f_twall(i_dim).eq.4) then
                appl_vel(i_dim)=0.
                write(*,*) "motion stopped at time step",turn_time(i_dim)
            end if
        end if

        if(f_twall(i_dim).eq.5) then
            appl_vel(i_dim)=va_spring_twall(i_dim)*cos(2*pi*i_time/turn_time(i_dim))

            if(mod(i_time,turn_time(i_dim)).eq.0) then
                write(*,*) "new period at time step",i_time
            end if
        end if


        if(f_twall(i_dim).eq.6) then
            if(i_time.lt.p_time+turn_time(i_dim)) then
                appl_vel(i_dim)=va_spring_twall(i_dim)*(1-cos(pi*(i_time-p_time)/turn_time(i_dim)))/2
            end if
            if(i_time.ge.(p_time+turn_time(i_dim))) then
                appl_vel(i_dim)=va_spring_twall(i_dim)
            end if
            if(i_time.lt.p_time) then
                appl_vel(i_dim)=0.0
            end if
            if(i_time.eq.p_time) write(*,*) "starting motion at timestep",i_time
            if(i_time.eq.(p_time+turn_time(i_dim)))  write(*,*) "final velocity reached at timestep",i_time
        end if

        if(f_twall(i_dim).eq.7) then
            if(i_time.lt.(p_time+turn_time(i_dim))) then
                appl_vel(i_dim)=va_spring_twall(i_dim)*cos(pi*(i_time-p_time)/turn_time(i_dim))
            end if
            if(i_time.ge.(p_time+turn_time(i_dim))) then
                appl_vel(i_dim)=-va_spring_twall(i_dim)
            end if
            if(i_time.lt.p_time) then
                appl_vel(i_dim)=va_spring_twall(i_dim)
            end if
            if(i_time.eq.p_time) write(*,*) "starting inversion at timestep",i_time
            if(i_time.eq.(p_time+turn_time(i_dim)))  write(*,*) "inversion performed at timestep",i_time
        end if


        if(f_twall(i_dim).eq.8) then
            if(i_time.lt.(p_time+turn_time(i_dim))) then
                appl_vel(i_dim)=va_spring_twall(i_dim)*(cos(pi*(i_time-p_time)/turn_time(i_dim))+1)/2
            end if
            if(i_time.ge.(p_time+turn_time(i_dim))) then
                appl_vel(i_dim)=0.0
            end if
            if(i_time.lt.p_time) then
                appl_vel(i_dim)=va_spring_twall(i_dim)
            end if
            if(i_time.eq.p_time) write(*,*) "stopping motion at timestep",i_time
            if(i_time.eq.(p_time+turn_time(i_dim)))  write(*,*) "stop performed at timestep",i_time
        end if

        if(f_twall(i_dim).eq.9) then
            appl_vel_bw(i_dim)=0.5*dt_2*(fbw(i_dim)/(0.5*n_wall)+va_spring_twall(i_dim))/mass(n_part)
            appl_vel_tw(i_dim)=0.5*dt_2*(ftw(i_dim)/(0.5*n_wall)-va_spring_twall(i_dim))/mass(n_part)

            !shiftfactor for external force

            shifttw(i_dim)=appl_vel_tw(i_dim) + appl_acc_tw(i_dim)
            shiftbw(i_dim)=appl_vel_bw(i_dim) + appl_acc_bw(i_dim)

        end if

    end do 
#   endif

! ----- Update positions with Velocity Verlet 
!print*,"aceleracion"
!print*,a
!print*,"velocidad"
!print*,v
!print*,"Fin positions before verlet"

    do i_part = 1 , n_mon_tot


#       if PINNED==1 
        if (i_part.eq.part_init_d+1) cycle   !excluding particle 3 and particle 4 from
        !thermostatisation, in case without spring
        if (i_part.eq.part_init_e+1) cycle
#       endif

!        r0(1,i_part) = r0(1,i_part) + dt*v(1,i_part) + 0.5*dt_2*a(1,i_part) 
!        r0(2,i_part) = r0(2,i_part) + dt*v(2,i_part) + 0.5*dt_2*a(2,i_part) 
!        r0(3,i_part) = r0(3,i_part) + dt*v(3,i_part) + 0.5*dt_2*a(3,i_part) 


        r0(:,i_part) = r0(:,i_part) + dt*v(:,i_part) + 0.5*dt_2*a(:,i_part) 


#       if STORE == 1 
! ----- Update unfolded coordinates
        r0_unfold(:,i_part) = r0_unfold(:,i_part) + dt*v(:,i_part) + 0.5*a(:,i_part)*dt_2
#       endif 

! ----- Update velocities  to the half of the interval 

!try        v_half(1,i_part) =  v(1,i_part) + 0.5*dt*a(1,i_part)
!try        v_half(2,i_part) =  v(2,i_part) + 0.5*dt*a(2,i_part)
!try        v_half(3,i_part) =  v(3,i_part) + 0.5*dt*a(3,i_part)

!        v(1,i_part) =  v(1,i_part) + 0.5*dt*a(1,i_part)
!        v(2,i_part) =  v(2,i_part) + 0.5*dt*a(2,i_part)
!        v(3,i_part) =  v(3,i_part) + 0.5*dt*a(3,i_part)

#ifdef BIDIMENSIONAL
        v(1,i_part) =  v(1,i_part) + 0.5*dt*a(1,i_part)
        v(2,i_part) = 0.0
        v(3,i_part) =  v(3,i_part) + 0.5*dt*a(3,i_part)
#else 
        v(1,i_part) =  v(1,i_part) + 0.5*dt*a(1,i_part)
        v(2,i_part) =  v(2,i_part) + 0.5*dt*a(2,i_part)
        v(3,i_part) =  v(3,i_part) + 0.5*dt*a(3,i_part)
#endif
    end do
!print*,r0
!print*,"Fin positions after verlet"

# if WALL == 1 /*explicit wall*/
 ! Update Wall positions 
 ! Claudio 2009: wall position works only en X coordinate. I use only f_t_wall(1)

        do i_wall=1,n_wall/2
             r_wall_equi(:,i_wall) = r_wall_equi(:,i_wall) + va_spring_twall(:)*dt
        end do
        do i_wall=n_wall/2+1,n_wall
             r_wall_equi(:,i_wall) = r_wall_equi(:,i_wall) - va_spring_twall(:)*dt
        end do
#endif 
#ifdef SHEARED
        do i_wall=n_wall/2+1,n_wall
            do i_dim =1,3
                if(f_twall(i_dim).lt.9) r_wall_equi(i_dim,i_wall) = r_wall_equi(i_dim,i_wall) - appl_vel(i_dim)*dt
                if(f_twall(i_dim).eq.9) r_wall_equi(i_dim,i_wall) = r_wall_equi(i_dim,i_wall) + shiftbw(i_dim)
            end do
            if(f_twall(i_dim).lt.9) then
                r0_twall(i_dim)= r0_twall(i_dim) + appl_vel(i_dim)*dt
                r0_bwall(i_dim)= r0_bwall(i_dim) - appl_vel(i_dim)*dt
            end if
            if(f_twall(i_dim).eq.9) then
                r0_twall(i_dim)= r0_twall(i_dim) + shifttw(i_dim)
                r0_bwall(i_dim)= r0_bwall(i_dim) + shiftbw(i_dim)
            end if
        end do
#endif




            do i_part = 1,n_mon_tot
#if SYMMETRY == 0
                do i_dim = 1,n_dim-1
# elif SYMMETRY == 1
                    do i_dim = 1,n_dim
#endif
                        if(r0(i_dim,i_part).gt.boundary(i_dim)) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                            mic_count(i_dim,i_part) = mic_count(i_dim,i_part) + 1
                        else if(r0(i_dim,i_part).lt.(0.)) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                            mic_count(i_dim,i_part) = mic_count(i_dim,i_part) - 1
                        end if
                    end do
            end do
#   if WALL == 1 
!  ----  PBC conditions for wall atoms 
        do i_wall = 1,n_wall
            do i_dim = 1,n_dim-1
                if(r_wall_equi(i_dim,i_wall).gt.(boundary(i_dim))) then
                    r_wall_equi(i_dim,i_wall)=r_wall_equi(i_dim,i_wall) - boundary(i_dim)
                    mic_count(i_dim,i_wall) = mic_count(i_dim,i_wall) + 1
                else if(r_wall_equi(i_dim,i_wall).lt.(0.)) then
                    r_wall_equi(i_dim,i_wall)=r_wall_equi(i_dim,i_wall) + boundary(i_dim)
                    mic_count(i_dim,i_wall) = mic_count(i_dim,i_wall) - 1
                end if
            end do
        end do
    
#endif
        !test
# if WALL == 1 /* explicit wall */        
! Claudio 2009: Positions are fixed to equilibrium sites if we have explicit
! walls
           do i_wall = 1,n_wall
               i_part=n_mon_tot+i_wall
               r0(:,i_part)=r_wall_equi(:,i_wall)
           end do
#   endif    
! --- Calling the geometric routine for reshaping the droplet
!
!DROP, WARN: reshape in each time step if thermalizing

#   if SYSTEM==1 /* droplet */
        if (i_time <= n_relax ) then
          call geom_constraint()
        end if
#   endif



!
!   *** DEBUGGING STUFF *** 
!

    if(debug) then
        print*, "random force= " ,sum ( sum (force_r(:,:) , dim=1) ,dim=1)
        print*, "dissip force= " ,sum ( sum (force_d(:,:) , dim=1) ,dim=1)
    end if
end subroutine verlet_positions
