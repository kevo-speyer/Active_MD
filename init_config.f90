     subroutine init_config
#include 'control_simulation.h' 
      use commons 
      use util 
      implicit none
      logical :: test_file
      logical,parameter :: debug=.false.
      
!---  xy plane

      x_max = x_space * dble(n_cell_w_x)
      y_max = y_space * dble(n_cell_w_y)
      surface = x_max*y_max

!---  z direction

      boundary(1) = x_max
      boundary(2) = y_max
      boundary(3) = z_space_wall
      inv_boundary(:) = 1./boundary(:)
      inv_V=inv_boundary(1)*inv_boundary(2)*inv_boundary(3)
      half_boundary(:) = boundary(:)/2.


#if SYSTEM == 1 /* Droplets */
!
! ***** Sizes and dimensions of droplet, vapour phase etc . *******
!

!
! --- Using the same units here define the dimensions of the droplet

      x_max_d = x_space * drop_cell_x
      y_max_d = y_space * drop_cell_y

!      
! Boundaries for the droplet box

      boundary_d(1) = x_max_d
      boundary_d(2) = y_max_d
      boundary_d(3) = z_max_d ! The droplet goes only the 70% of the position of the top wall
!      
! Dimensions of droplet box in x
!

      x_box_min = boundary(1)/2. - boundary_d(1)/2.  
      x_box_max = boundary(1)/2. + boundary_d(1)/2.

! ***** Droplet stuff  ****************************

!
! Get parameters of the droplet and write them out
!
!---------------------------------------------
!NOTE:*  z_min_d should be the estimated high of the brush polymer layers.
!       It is only used for the estimation of the density and volume of the droplet.
!      Around this value is expected the drop/brush interface
!     *  z_max_d is the real high of the droplet. During thermalization the
!      z of droplets chains will be confined up to this value.
!---------------------------------------------

   ! Estimated droplet density     

   vol_drop = x_max_d*y_max_d*(z_max_d-z_min_d)
   dens_drop =  dble(n_mon_d*n_chain_d)/vol_drop   

   print '(/a/)'," **** Droplet data ****"
   print '(/a,3(f8.5,1x))',"    Droplet dims (z estimated)= ",x_max_d,y_max_d,(z_max_d-z_min_d) 
   print '(a,2f8.5)',"    From system_input z_min_d,z_max_d= ",z_min_d,z_max_d
   print '(/a,f15.9)',"   Droplet vol.=", vol_drop
   print '(a,f15.9/)',"   Estimated droplet density= ",dens_drop 

#endif /*system=1; droplets*/

! -----  Size of binning boxes and allocation of variables     

        call make_binning_boxes


! Checks existence of the old config file

      inquire(file="conf_old",exist=test_file)

      if(.not.test_file) then  ! --------        IF conf_old DOES NOT EXIST 

          print '(/a/)'," *   'conf_old' file not found. Using default configuration ..."

! *** Generate default configuration

          call conf_default()


#if PINNED==1|| PINNED==2
      r0(3,part_init_d+1)= boundary(3)/2
      r0(3,part_init_e+1)= boundary(3)/2
      r0(2,part_init_e+1)= r0(part_init_d+1, 2)
      r0(1,part_init_e+1)= r0(part_init_d+1, 1)+ dist_pinned
      r_start(1,:)=r0(:,part_init_d+1)
      r_start(2,:)=r0(:,part_init_e+1)

#endif

      else  ! ---------------------     IF conf_old  EXIST -------------------------

!
! ----- Read conf_old
!

          open(10,file="conf_old",status="old")
          write(6,'(/a)',advance='no')  "  *   Reading  conf_old..."
! ----     Check whether number of monomers right
          read(10,*) !info !(antes )
#ifdef STARS

          read(10,'(2i10)') i_mon,i_chain
          read(10,'(3i10)') i_mon_d,i_chain_d
          read(10,'(4i10)') i_mon_arm,i_arm,i_star
          if((i_mon.ne.n_mon).or.(i_chain.ne.n_chain)) then
              write(*,'(a)') "   Polymer set-ups not compatible"
              if((i_mon*i_chain).ne.(n_mon*n_chain)) then
                  print *, "    Warning: n_mon and/or n_chain of mfa_input"
                  print *, "and system input don't coincide !"
                  print *," This may be a problem if using confs from conf_old"
                  !               print*,"Stopping here ..."
                  stop
              end if
              if((i_mon_d.ne.n_mon_d).or.(i_chain_d.ne.n_chain_d)) then
                  writE(*,'(a)') "Solvent set-ups not compatible"
                  stop
              end if
              if((i_mon_arm.ne.n_mon_arm).or.(i_arm.ne.n_arms).or.(i_star.ne.n_stars)) then
                  writE(*,'(a)') "Stars set-ups not compatible"
                  stop
              end if

          end if
#else
      read(10,'(2i10)') i_mon,i_chain
          if((i_mon.ne.n_mon).or.(i_chain.ne.n_chain)) then
              write(*,'(a)') "  * Polymer set-ups not compatible"
              if((i_mon*i_chain).ne.(n_mon*n_chain)) then
                  print *, "  *   Warning: n_mon and/or n_chain of mfa_input and system input don't coincide !"
                  print *,"   *   This may be a problem if using confs from conf_old"
              end if
         end if
#endif
          read(10,120) i_cell_w_x,i_cell_w_y !we read in any case just for compatibility 

#       if WALL == 1
!          if(f_explicit_wall) then
!----  check whether number of wall unit cells right
              if((i_cell_w_x.ne.n_cell_w_x).or.(i_cell_w_y.ne.n_cell_w_y)) then
                  write(*,'(a)') "Set-up of wall not compatible"
                  if((i_cell_w_x*i_cell_w_y).ne.(n_cell_w_x*n_cell_w_y)) stop
              end if
!          end if       
#       endif

          read(10,110) j_order
          j_order = min(j_order,n_order)

!----  get configuration and its derivatives

          do i_part = 1,n_part
              read(10,113) a_type(i_part),r0(1,i_part),r0(2,i_part),r0(3,i_part)
          end do

          do i_part = 1,n_part
              read(10,103) (v(i_dim,i_part),i_dim=1,n_dim)
          end do
          do i_part = 1,n_part
              read(10,103) (a(i_dim,i_part),i_dim=1,n_dim)
          end do

! Note: the program writes half of acceleration and here reads half of accelerations
! coherently
! ------
!          if(f_explicit_wall) then !if no explicit wall doesn't read wall atoms
#if     WALL == 1 /* explicit wall */
!----  equilibrium sites of wall atoms
              do i_wall = 1,n_wall
                  read(10,103) (r_wall_equi(i_dim,i_wall),i_dim=1,n_dim)
                  i_part = i_wall + n_mon_tot
                  r0(:,i_part) = r_wall_equi(:,i_wall)
              end do
!----  displacement of top wall and its derivatives
              read(10,103) (r0_twall(i_dim),i_dim=1,n_dim)
              do i_order = 1,j_order
                  r_dummy = dt**i_order
                  read(10,103) (rx_twall(i_dim,i_order),i_dim=1,n_dim)
                  do i_dim = 1,n_dim
                      rx_twall(i_dim,i_order) = rx_twall(i_dim,i_order)*r_dummy
                  end do
              end do


              read(10,130) (pbc_twall(i_dim),i_dim=1,n_dim)
!----  displacement of bottom wall
              read(10,103) (r0_bwall(i_dim),i_dim=1,n_dim)
              read(10,130) (pbc_bwall(i_dim),i_dim=1,n_dim)
              read(10,103) (r0_spring_twall(i_dim),i_dim=1,n_dim)
              read(10,130) (s_force_grad(i_dim),i_dim=1,n_dim) 
#       endif
!          end if  ! explicit wall atoms     

          print '(a/)', "        Done !"

          close(10) ! End configuration file reading
      
#if PINNED==1|| PINNED==2
      r_start(1,:)=r0(:,part_init_d+1)
      r_start(2,:)=r0(:,part_init_e+1)
#endif
      end if !  ----------------------- END IF conf_old ------------------

! HERE: either conf_old was read or first conf was generated


#if SYSTEM == 2  || SYSTEM == 3
!
! Parameters  for  Coulomb interactions
!
    q(:)=0.
    do i_part = 1 ,n_mon_tot  !!! walls are excluded, they will not be charged
                              !!! coulombic interaction should not be calculated
                              !!!with them, only LJ
       q(i_part) = charge(a_type(i_part))
    end do
! --- Initialize charges for all the particles 

    sqrt_alpha = sqrt(alpha_c)
    inv_sqrt_pi = 1./pi


! Self interaction calculation for the Ewald sum 
    
      v_coul_self = 0. 
      do i_part = 1 , n_mon_tot 
      v_coul_self = v_coul_self - q(i_part)**2
      end do
      v_coul_self = v_coul_self*sqrt_alpha/sqrt(pi)


! Initialize Ewald variables in the k-space 

        call ewald_k(0)
!     lz_half = ( boundary(3)+z_slab*boundary(3) )/2.
     lz_half =  boundary(3)/2.
#endif /*system=2: coulomb*/


! ---- First binning and check fluid-fluid interactions

#       if BIN_TYPE == 0 
          call binning()
#       elif BIN_TYPE == 1         
          call my_binning()
#       endif        
  
          call fluid_fluid()
          print '(/a,f16.5/)',"    * V_fluid_fluid for the first configuration = ",v_fluid_fluid
          force(:,:)= 0.
          call check_fluid_fluid()

!
! Force switch on even with conf_old if s_time = 1        
!
#       ifdef FORCE_SWITCH_ON
          if(s_time == 1) then
              print '(/a/)',"  *  s_time=1. The simulation will use force switch_on " 
              call  init_force_switch_on()
          else
!             
!----- Make no ramp for r_2_min (force switch on not used)
!
              print '(/a)',"  *  s_time not 1 and conf_old exist. The simulation will NOT  use force switch_on " 
              print '(a/)',"     (Real force from start up)"

              r_2_min_init = 0.
              r_2_min_time = 0.
              r_2_min = r_2_min_init !Actually, it was already deifined in
                                     !fluid_fluid 
                                     !if run is a follow up, no switching on
          end if ! s_time = 1
#       endif
!

!
!---  reinitialize velocity of top wall if necessary
!
#       if WALL == 1 
          do i_dim = 1,n_dim
              if(f_twall(i_dim).eq.1) then
                  rx_twall(i_dim,1) = va_spring_twall(i_dim)*dt
                  do i_order = 2,n_order
                      rx_twall(i_dim,i_order) = 0.
                  end do
              end if
          end do
#       else    
          rx_twall(:,:) =0.0
#       endif

!
!
! -----  store masses and their inverses in array

      do i_part = 1,n_part
          mass(i_part)     = mass_type( a_type(i_part) )
          inv_mass(i_part) = 1./mass_type( a_type(i_part) )
      end do

! this is 0 if not explicit wall
      mass_twall = 4*n_wall*mass_type(n_type)
!
! First CM velocity   ! WARNING: here is assummed that the heads already have v=0

      vec_dummy(:) = 0.
      r_dummy = 0.
      do i_part = 1 , n_mon_tot 
          vec_dummy(:) = vec_dummy(:) + mass(i_part)*v(:,i_part)
          r_dummy = r_dummy + mass(i_part)
      end do
      vec_dummy(:) = vec_dummy(:)/r_dummy
      print '(/a,3f16.5/)',"   - First Vcm=", vec_dummy(:) 


#       if SYMMETRY == 1 
        print '(/a/)', "  *   Substrating Vcm so that initial tot. P=0 " 
        do i_part=1,n_mon_tot
            v(1,i_part) = v(1,i_part) - vec_dummy(1) 
            v(2,i_part) = v(2,i_part) - vec_dummy(2) 
            v(3,i_part) = v(3,i_part) - vec_dummy(3) 
        end do
#       endif        

      vec_dummy(:) = 0.
      r_dummy = 0.
      do i_part = 1 , n_mon_tot 
          vec_dummy(:) = vec_dummy(:) + mass(i_part)*v(:,i_part)
          r_dummy = r_dummy + mass(i_part)
      end do
      vec_dummy(:) = vec_dummy(:)/r_dummy
      print '(/a,3f16.5/)',"   - First Vcm(substracting Vcm)=", vec_dummy(:) 

#ifdef POISEUILLE
! The readed value is a gravity g=f/m. The actual  force must be scaled with the
! mass of each particle 

! Allocation of the external force vector

       allocate(ext_force(3,n_mon_tot)) 

! Define the value
        ext_force(:,:) = 0.
        do i_part = 1,n_mon_tot
            if (a_type(i_part) /= 1 ) then ! ext force is not applied in the brush head 
                ext_force(:,i_part) = mass(i_part)*const_force(:)
            end if
        end do

#ifdef PARTICLE_4
       print '(/a,3f15.10/)' ," Const. External 'gravity' =", const_force(:) !  [POSEUILLE VERSION]
       print '(a,4f12.8)' ," Const. External Fx =", 0.00,const_force(1)*(mass_type(2:)) !  [POSEUILLE VERSION]
       print '(a,4f12.8)' ," Const. External Fy =", 0.00,const_force(2)*(mass_type(2:)) !  [POSEUILLE VERSION]
       print '(a,4f12.8/)' ," Const. External Fz =", 0.00,const_force(3)*(mass_type(2:)) !  [POSEUILLE VERSION]
#else
       print '(/a,3f15.10/)' ," Const. External 'gravity' =", const_force(1:3) !  [POSEUILLE VERSION]
       print '(a,4f12.8)' ," Const. External Fx =", 0.00,const_force(1)*(mass_type(2:3)) !  [POSEUILLE VERSION]
       print '(a,4f12.8)' ," Const. External Fy =", 0.00,const_force(2)*(mass_type(2:3)) !  [POSEUILLE VERSION]
       print '(a,4f12.8/)' ," Const. External Fz =", 0.00,const_force(3)*(mass_type(2:3)) !  [POSEUILLE VERSION]
#endif
#endif

#ifdef SHEARED

        write(*,*)
        write(*,*) "Shear protocol:"

if (f_twall(1).eq.0) then
        write(*,*) "x-direction: apply constant velocity (instanteneous)",va_spring_twall(1),"after time step",turn_time(1)
end if


if (f_twall(1).eq.1) then
        write(*,*) "x-direction: constant velocity",va_spring_twall(1)
end if

if (f_twall(1).eq.2) then
        write(*,*) "x-direction: invert shear direction after step (instanteneous)",turn_time(1) 
end if

if (f_twall(1).eq.3) then
        write(*,*) "x-direction: oscillatory motion (instantaneous) with period",turn_time(1) 
end if

if (f_twall(1).eq.4) then
        write(*,*) "x-direction: stop motion (instanteneous) at time step",turn_time(1)
end if

if (f_twall(1).eq.5) then
        write(*,*) "x-direction: oscillatory motion (sinosoidal) with period",turn_time(1) 
end if

if (f_twall(1).eq.6) then
        write(*,*) "x-direction: apply constant velocity after",p_time,"steps within",turn_time(1),"steps" 
end if

if (f_twall(1).eq.7) then
        write(*,*) "x-direction: invert shear direction after",p_time,"steps within",turn_time(1),"steps" 
end if

if (f_twall(1).eq.8) then
        write(*,*) "x-direction: stop motion after",p_time,"steps within",turn_time(1),"steps" 
end if

if (f_twall(1).eq.9) then
        write(*,*) "x-direction: apply external force per wall atom",va_spring_twall(1);
end if



if (f_twall(2).eq.0) then
        write(*,*) "y-direction: apply constant velocity (instanteneous)",va_spring_twall(2),"after time step",turn_time(2)
end if


if (f_twall(2).eq.1) then
        write(*,*) "y-direction: constant velocity",va_spring_twall(2)
end if

if (f_twall(2).eq.2) then
        write(*,*) "y-direction: invert shear direction (instanteneous) after step",turn_time(2) 
end if

if (f_twall(2).eq.3) then
        write(*,*) "y-direction: oscillatory motion (instantaneous) with period",turn_time(2) 
end if

if (f_twall(2).eq.4) then
        write(*,*) "y-direction: stop motion (instanteneous) at time step",turn_time(2)
end if

if (f_twall(2).eq.5) then
        write(*,*) "y-direction: oscillatory motion (sinosoidal) with period",turn_time(2) 
end if

if (f_twall(2).eq.6) then
        write(*,*) "y-direction: apply constant velocity after",p_time,"steps within",turn_time(2),"steps" 
end if

if (f_twall(2).eq.7) then
        write(*,*) "y-direction: invert shear direction after",p_time,"steps within",turn_time(2),"steps" 
end if

if (f_twall(2).eq.8) then
        write(*,*) "y-direction: stop motion after",p_time,"steps within",turn_time(2),"steps" 
end if

if (f_twall(2).eq.9) then
        write(*,*) "y-direction: apply external force per wall atom",va_spring_twall(2);
end if


if (f_twall(3).eq.0) then
        write(*,*) "z-direction: apply constant velocity (instanteneous)",va_spring_twall(3),"after time step",turn_time(3)
end if

if (f_twall(3).eq.1) then
        write(*,*) "z-direction: constant velocity",va_spring_twall(3)
end if

if (f_twall(3).eq.2) then
        write(*,*) "z-direction: invert shear direction (instanteneous) after step",turn_time(3)
end if

if (f_twall(3).eq.3) then
        write(*,*) "z-direction: oscillatory motion (instantaneous) with period",turn_time(3) 
end if

if (f_twall(3).eq.4) then
        write(*,*) "z-direction: stop motion (instanteneous) at time step",turn_time(3)
end if

if (f_twall(3).eq.5) then
        write(*,*) "z-direction: oscillatory motion (sinosoidal) with period",turn_time(3) 
end if

if (f_twall(3).eq.6) then
        write(*,*) "z-direction: apply constant velocity after",p_time,"steps within",turn_time(3),"steps" 
end if

if (f_twall(3).eq.7) then
        write(*,*) "z-direction: invert shear direction after",p_time,"steps within",turn_time(3),"steps" 
end if

if (f_twall(3).eq.8) then
        write(*,*) "z-direction: stop motion after",p_time,"steps within",turn_time(3),"steps" 
end if

if (f_twall(3).eq.9) then
        write(*,*) "z-direction: apply external force per wall atom",va_spring_twall(3);
end if
write(*,*)
#endif
      if(debug) then
          vec_dummy(:) = 0.
          do i_part=1,n_part
              if(r0(3,i_part)<0.5 .or. r0(3,i_part)>(z_space_wall-0.5))  &
              print*,"[init_config] Particle too close to the wall",r0(3,i_part) 
              write(100,'(9f10.5)') r0(:,i_part),v(:,i_part),a(:,i_part)
              vec_dummy(:) = mass(i_part)*v(:,i_part)
          end do
!stop
          print'(a30,3f15.5)',"   - First P=",vec_dummy(:)
      end if

#ifdef STARS   

      103 format(3f25.7)
#else
      103 format(3f15.7)
#endif
      110 format(i3)
      113 format(i3,3f13.7)
      120 format(2i10)
      130 format(3i13)
      202 format(3f15.6,a)
      210 format(a,i3)
      220 format(3i15,a)
 end subroutine init_config
