subroutine conf_default()
#include 'control_simulation.h'
! Generates the default configuration for
! - wall atoms
! - brushes or adsorbed polymers
! - explicit solvent monomers (probably not anymore)
! - the first droplet generation or melt
! - Initialize here the force "switch-on" 

      use commons 
      use util 
      use ziggurat
!      use functions
      implicit none

      logical, parameter :: debug=.false.
!      real (kind=8) :: z_head 
! distance of brushes head from the wall 

#      if     WALL == 1 /* explicit wall */
          z_head =  2.**(1./6.)    !original value. It is the same value !!!
#       else
          z_head = 1.2*sigma(1,1) ! 1.2sigma of brushes heads and body 
#       endif
!      end if

!      do i_dim = 1,n_dim
!          do i_part = 1,n_part
              r0(:,:) = 0.
!          end do
          r0_twall(:) = 0.
          s_force_grad(:) = 1
!      end do
      r0_twall(3) = z_space_wall

! ------------ LABELLING OF PARTICLES 

!---  define which particle belongs to which species
!! debug
#if SYSTEM == 0 || SYSTEM == 1  || SYSTEM == 4
!
! ---  Brush build-up 
!
      do i_part = 1,part_init_d
       if(mod(i_part-1,n_mon).eq.0) then
       a_type(i_part) = 1   !first monomer of the chain
       else
       a_type(i_part) = 2 ! the rest of the chain
       endif
      end do

#ifndef PARTICLE_4
! Melt or droplet 
      do i_part = part_init_d+1,n_mon_tot
       a_type(i_part) = 3             
      end do
#else
! Melt or droplet
      do i_part = part_init_d+1,part_init_e
       a_type(i_part) = 3               
      end do
! Particle 4
      do i_part = part_init_e+1,n_mon_tot
       a_type(i_part) = 4              
      end do
#endif /* Particle 4 */

#elif SYSTEM == 2  /* Charged system */

! Ions with charge +
      do i_part = 1,part_init_d
       a_type(i_part) = 2   !first monomer of the chain
      end do
! Ions with charge -
      do i_part = part_init_d+1,part_init_e
       a_type(i_part) = 3             
      end do
! Particle 4: Now charged polymers
      do i_part = part_init_e+1,n_mon_tot
       a_type(i_part) = 4              
      end do

#elif SYSTEM == 3

!
! ---  Brush build-up 
!
      do i_part = 1,part_init_d
       if(mod(i_part-1,n_mon).eq.0) then
       a_type(i_part) = 1   !first monomer of the chain
       else
       a_type(i_part) = 2 ! the rest of the chain
       endif
      end do

! Ions, polymer melt
      do i_part = part_init_d+1,part_init_e
       a_type(i_part) = 3             
      end do
! Ions, polymer melt
#       ifdef PARTICLE_4
             do i_part = part_init_e+1,n_mon_tot
              a_type(i_part) = 4              
             end do
#       endif

#endif


!    if(f_explicit_wall) then !all this if we consider explicit wall atoms
#   if WALL == 1 
        call gen_wall()
#  endif
!    end if ! explicit wall


! ----Build up FIRST CONFIGURATION

! Note: the algorithm doesn't work for only one brush

! ---- Build up brushes 

#if SYMMETRY == 0 
#   if SYSTEM ==0   /* Brush un top and bottom walls: channel geometry */
#       if BRUSH_TYPE ==1  /* Avoids overlp between grafting beads*/
            call gen_brush(3) 
#       elif BRUSH_TYPE ==0              /* Allows overlp between grafting beads*/
            call gen_brush(1)
#       elif BRUSH_TYPE ==2              /* Ordered Brush*/
            call gen_brush(5)   
#       elif BRUSH_TYPE ==3              /* Ordered Brush, polymers aligned*/
            call gen_brush(7)
#       endif
#   elif SYSTEM==1  /* Brush only in bottom wall: droplet */
#       if BRUSH_TYPE ==1  /* Avoids overlp between grafting beads*/
            call gen_brush(4)
#       elif BRUSH_TYPE ==0              /* Allows overlp between grafting beads*/
            call gen_brush(2)
#       elif BRUSH_TYPE ==2              /* Ordered Brush*/
            call gen_brush(6)
#       elif BRUSH_TYPE ==3              /* Ordered Brush, polymers aligned*/
            call gen_brush(8) !not working!!!
#       endif
#   endif /* a system with brush chains */
#endif  /*SYMMETRY*/


! ---- Call the drop/melt generating routine.
! 
#if SYSTEM == 0  || SYSTEM == 4 
      call gen_droplet(1)
#       ifdef PARTICLE_4
!        Defines the positions of particle 4
             call gen_droplet(2)
#       endif
#endif


#if SYSTEM == 1
! Generate the droplet, from a cubic shape
      call gen_droplet(4)
#endif

#       if SYSTEM == 2
!        Generate particle 2 coordinates when there are no brushes 
             call gen_droplet(3)
!         Generates particle 3: ions  or polymers 
             call gen_droplet(1)
#       ifdef PARTICLE_4
!        Generates particle 4: ions or polymers
             call gen_droplet(2)
#       endif
#       endif
#if SYSTEM == 3
! Brush in lower wall 
           call gen_brush(1)  ! particle types 1 AND 2 

!        Generate particle 2 coordinates when there are no brushes 
!             call gen_droplet(3)
!         Generates particle 3: ions  or polymers 
             call gen_droplet(1)
#       ifdef PARTICLE_4
!        Generates particle 4: ions or polymers
             call gen_droplet(2)
#       endif


#endif
#ifdef STARS
        call gen_droplet(5)     !generates stars, the types are also given there
#endif

! debug
! call write_conf(2,r0,10)

! ***********************  Setup 'force switch-on':  the configuration is random
!
!
#   ifdef FORCE_SWITCH_ON
      call  init_force_switch_on()
#   endif
!
!---  Initialize accelerations
!
#   if BIN_TYPE == 0 
       call binning()
#   elif BIN_TYPE == 1 
       call my_binning()
#   endif

       force(:,:) = 0.
       call fluid_fluid()

#if SYMMETRY == 0
#       if WALL !=1 /* not explicit wall*/
         call fluid_wall(wall_flag) ! 2= wall potential 9-3 (without atoms)
#       endif                      
!try without                        call wall_wall(wall_flag)
#   if WALL == 1 /* explicit wall */
!try without                        call intra_wall
#   endif
#endif


        call intra_molec()

 print '(//a,f16.5)',"   - First v_wall_wall= ",v_wall_wall 
 print '(a,f16.5)',  "   - First v_fluid_wall= ",v_fluid_wall 
 print '(a,f16.5)',  "   - First v_fluid_fluid= ",v_fluid_fluid 
 print '(a,f16.5//)',"   - First Mean square force per particle= ", sum( sum( force**2 , dim=1),dim=1)/real(n_mon_tot) 
 
!---  Convert forces into accelerations
!
      mass_twall = 4*n_wall*mass_type(n_type)
      mass_twall = mass_twall/2
      !old: n_loop
      do i_part = 1,n_part
          do i_dim= 1,n_dim
              a(i_dim,i_part) = force(i_dim,i_part)/mass_type(a_type(i_part)) 
          end do
      end do
!
! ---- Generate a first set of initial velocities 
! 

#if VEL_INIT == 0         
        print *,"  * Setting initial velocities to ZERO "

        do i_part = 1 , n_mon_tot
            v(1,i_part)  = 0.0
            v(2,i_part)  = 0.0 
            v(3,i_part)  = 0.0
        end do
 
#elif VEL_INIT == 1         
        print *,"  *  Generating Maxwell-Boltzmann distributed initial velocity "

        do i_part = 1 , n_mon_tot
            v(1,i_part)  =  sqrt(temp/mass_type(a_type(i_part)))*rnor()
            v(2,i_part)  =  sqrt(temp/mass_type(a_type(i_part)))*rnor()
            v(3,i_part)  =  sqrt(temp/mass_type(a_type(i_part)))*rnor()
        end do 
#endif

#ifdef BIDIMENSIONAL
        v(2,:) = 0.0
#endif

#if WALL == 1 /* explicit wall*/
           rx_twall(:,2) = force_twall(:)/mass_twall *0.5*dt_2
!       else
#else /*every other wall */
           rx_twall(:,2)=0. ! if not exlicit wall atoms, set accel wall = 0
#endif
       
! Set to  forces in brushes if applies

! If not explicit wall: we cancel the force over the brushes heads even for the first configutrations
! we do that overriding the previous calculation

!
#if SYMMETRY == 0 
#   if SYSTEM == 0 || SYSTEM == 1  
#       if WALL == 2 || WALL == 3 || WALL ==4 
        if(n_chain>0) then ! the wall force compensates always the force from the system over the brushes
            force(:,1:n_chain*n_mon:n_mon) = 0. !force
                a(:,1:n_chain*n_mon:n_mon) = 0.  ! accel
                v(:,1:n_chain*n_mon:n_mon) = 0.  ! velocities 
        end if
#       endif
#   endif
#endif

    if(debug) then
#if SYMMETRY == 0         
              do i_part=1,n_mon_tot
                  if(r0(3,i_part)>boundary(3) - 0.5 .or.r0(3,i_part)<0.5) then
                      print *, "[conf_default] Wall_too_close, z =",r0(3,i_part)
                  end if
              end do
#endif              
    end if

      !call write_conf(1,r0,10)  ; stop

      end  subroutine conf_default 
