subroutine conf_default()
 

 




 



 
 



 

 



 

 












 




                        




 




 


                     





                     








                  
                  
                  
 

    


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


          z_head = 1.2*sigma(1,1) ! 1.2sigma of brushes heads and body

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


! Melt or droplet
      do i_part = part_init_d+1,n_mon_tot
       a_type(i_part) = 3             
      end do





!    if(f_explicit_wall) then !all this if we consider explicit wall atoms

!    end if ! explicit wall


! ----Build up FIRST CONFIGURATION

! Note: the algorithm doesn't work for only one brush

! ---- Build up brushes




            call gen_brush(7)   





! ---- Call the drop/melt generating routine.
!

      call gen_droplet(1)










! debug
! call write_conf(2,r0,10)

! ***********************  Setup 'force switch-on':  the configuration is random
!
!

      call  init_force_switch_on()

!
!---  Initialize accelerations
!

       call binning()


       force(:,:) = 0.
       call fluid_fluid()



         call fluid_wall(wall_flag) ! 2= wall potential 9-3 (without atoms)

!try without                        call wall_wall(wall_flag)




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
! ---- Generate initial velocities
!
         

        print *,"  *  Setting initial velocities to ZERO "


        do i_part = 1 , n_mon_tot

            v(1,i_part)  =  0.
            v(2,i_part)  =  0.
            v(3,i_part)  =  0.

        end do

 

        v(2,:)= 0.

        



           rx_twall(:,2)=0. ! if not exlicit wall atoms, set accel wall = 0

       
! Set to  forces in brushes if applies

! If not explicit wall: we cancel the force over the brushes heads even for the first configutrations
! we do that overriding the previous calculation

!



        if(n_chain>0) then ! the wall force compensates always the force from the system over the brushes
            force(:,1:n_chain*n_mon:n_mon) = 0. !force
                a(:,1:n_chain*n_mon:n_mon) = 0.  ! accel
                v(:,1:n_chain*n_mon:n_mon) = 0.  ! velocities
        end if




    if(debug) then

              do i_part=1,n_mon_tot
                  if(r0(3,i_part)>boundary(3) - 0.5 .or.r0(3,i_part)<0.5) then
                      print *, "[conf_default] Wall_too_close, z =",r0(3,i_part)
                  end if
              end do

    end if

!call write_conf(1,r0,10)  ; stop

      end  subroutine conf_default 
