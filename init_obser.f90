      subroutine init_obser
      use commons 
      use functions 

#include 'control_simulation.h'
! OBSOLETE #ifndef PROFILES
! OBSOLETE ! conflicts with vel_prof from functions.f90 
! OBSOLETE       use util , only: velocity_prof
! OBSOLETE #else
#ifdef PROFILES
      use functions
#endif      
      implicit none
      integer :: i_chain_m
      character (len=4) , parameter , dimension(8) :: nam1= (/"Cl  ","O   ","He  ","N   ","S   ","Al   ","H   ","Cu   "/)
      logical :: ufile 
      real (kind=8) :: rgx
!

      count_obs = 0
!
      v_wall_wall_1 = 0.
      v_wall_wall_2 = 0.
      t_wall_1 = 0.
      t_wall_2 = 0.
      e_wall_wall_1 = 0.
      e_wall_wall_2 = 0.
!
      v_fluid_fluid_1 = 0.
      v_fluid_fluid_2 = 0.
      t_fluid_1 = 0.
      t_fluid_2 = 0.
      e_fluid_fluid_1 = 0.
      e_fluid_fluid_2 = 0.
!
      v_fluid_wall_1 = 0.
      v_fluid_wall_2 = 0.
!
      v_total_1 = 0.
      v_total_2 = 0.
      t_total_1 = 0.
      t_total_2 = 0.
      e_total_1 = 0.
      e_total_2 = 0.
!clau: intra
      v_tot_intra_1 =   0.      
      v_tot_intra_2 =   0.
      v_tot_intra_d_1 = 0.
      v_tot_intra_d_2 = 0.

! Total momentum of the system 

      tot_P(:)   = 0.
      tot_P_2(:) = 0.


      time_ave_count = 0
      time_ave_count_2 = 0

!! obso      if(f_logscale.eq.0) then
!! obso          n_time_ave = 200 !if not log scale counts first 200 steps before calculating things in obser
!! obso      end if
!! obso      if(f_logscale.eq.1) then
!! obso       n_time_ave = 0
!! obso      end if

      do i_dim = 1,n_dim
       r0_twall_1(i_dim) = 0.
       f0_spring_1(i_dim) = 0.
      end do
!
      n_fcell_max = 0
      n_wcell_max = 0
!
      n_neigh_fl_max = 0
      n_neigh_fw_max = 0
      n_neigh_ww_max = 0
!
      do i_dim = 1,n_dim
       r0_last_pinned(i_dim) = r0_twall(i_dim)
      end do
!
      z_1 = 0.
      z_2 = 0.
! Radius of Gyration mean variables
        r_g2_mean(:,:) = 0.0
        r_par_per2(:) = 0.0
! End to end radius:
        re_2_mean(:,:) = 0.0
! Force on heads variable 
        f_on_heads(:,:) = 0.0 


! Initialization of store_config routine
 
 call store_config(1)
 
 ! *** For diffusive quantities, drop the initial reference configuration 
 !      of the run (after a first equilibration run)
    
      if( s_time == 1 ) then
print '(/a/)',"  *  Writing the very first generated configuration to conf0.xyz file" 
        open(unit=23,file='conf0.xyz',status='unknown')
          write(23,*) n_part !ORI n_mon_tot
          write(23,*) "    First configuration" 
!
          do i_part = 1 , n_part
          write(23,'(a4,3f15.6)') nam1(a_type(i_part)),r0(:,i_part) 
          end do
          write(23,*) 

        close(23)
      end if

#ifdef PROFILES
      print  '(/a/)',"  *  Calculating profiles during the run" 
!
!  --- Parameters
!
      dim_prof_dens = 200 
      inv_n_mon = 1./dble(n_mon)
      inv_n_chain = 1./dble(n_chain)
      
#if SYSTEM == 0 || SYSTEM == 2 || SYSTEM == 3 /* channel or charges*/
!  --- Allocation       
!
#ifndef PARTICLE_4
        allocate ( v_prof(dim_prof_dens,3,2),v_prof_2(dim_prof_dens,3,2),histo_out(dim_prof_dens,2) )
#else
        allocate ( v_prof(dim_prof_dens,3,3),v_prof_2(dim_prof_dens,3,3),histo_out(dim_prof_dens,3) )
#endif
!
!  --- Init
!
!       call histo(0,dim_prof_dens,n_mon,n_chain,n_mon_d,n_chain_d,surface,boundary(3)) !,histo_out) 

        call histo(0,dim_prof_dens)

       zz_min = 0.

       call  velocity_prof(0,dim_prof_dens,zz_min,boundary(3)) !,v_prof(:,:,:),v_prof_2(:,:,:))

#endif /*channel or charges*/

#if SYSTEM==1 /*droplet*/

! Init droplet measuring variables

        vcm_d(:) = 0.   ! droplet CM velocity
        vcm_d2(:)= 0.
 
! Files for droplet CM's position and velocity 
! common var drop_cm = droplet CM position 

! Find a good value for the droplet CM
!  Also defined the first value for N_CM (which third of the box will contain the CM)

       call rebuild_drop2(0)
       call rebuild_drop2(1)

       guessed_drop_cm(:) = drop_cm(:) 
       print '(a,i/)' ,"        DROP CM should be in the third number = ",n_cm 

!       call write_conf(2,r0,10)
!       stop

! ******************************* IMPORTANT ********************************************

! Profile and initialization from mide_drop_force

! Init. size variables for profiles'

        Lx_o_2(1) = boundary(1) ! For translation of the droplet to comoving frame 
        Lx_o_2(2:3) = 0.

        binning_box_length = 1. ! approx. dimension of the binning box in units of sigma
        
        n_box(2:3) = int (boundary(2:3)/binning_box_length) ! number of boxes in each direction
        n_box(1)=   int (2.*boundary(1)/binning_box_length) ! 
        
        n_box_xz(:) = n_box(:) ; n_box_xz(2) = 1

! Note:  In x dir, we add  boxes needed because of the unfolding. We double the number of boxes  

! Range covered by the binning boxes of the histograms [DEFINED. AVOID CHANGING THEM]. 
! They are function of n_box

        r_box_xz(1) = 2.*boundary(1)/dble(n_box(1)) ! For proper copy-pasting, the range in x is twice the box size
        r_box_xz(2) = boundary(2)
        r_box_xz(3) = boundary(3)/dble(n_box(3))

! Allocate vectors 

      allocate ( dens_xz(n_box(1),1,n_box(3)) , dens_xz_2(n_box(1),1,n_box(3)) )
      allocate (  dens_xz_b(n_box(1),1,n_box(3)) ,  dens_xz_b_2(n_box(1),1,n_box(3)) )
      allocate (  r0_unfold(3,n_part), r_cms(n_chain+n_chain_d,3) )
      allocate ( vel_xz(n_box(1),n_box(3),3), vel_xz_b(n_box(1),n_box(3),3)  )
      allocate ( vel_xz2(n_box(1),n_box(3),3), vel_xz_b2(n_box(1),n_box(3),3)  )

!! NOTE YET #ifdef FORCE_PROFILE      
!! NOTE YET       allocate ( force_prof(n_box(1),3) )
!! NOTE YET #endif

! Init. routines (all in util.f90 )

      call dens_prof(0,dens_xz(:,:,:),dens_xz_b(:,:,:),dens_xz_2(:,:,:),dens_xz_b_2(:,:,:),n_box_xz, &
      r_box_xz) ! end 
        
      call vel_prof_2d(0,vel_xz,vel_xz2,vel_xz_b,vel_xz2,n_box_xz,r_box_xz,dens_xz) ! 


        print '(/a,3f10.5)', "First drop CM = " ,  guessed_drop_cm(:)
        print '(a/)', "  **** CHECK IF IT IS NOT TOO BAD ****** " 

!  Initialize 

      call calc_cms(1) 

      call rebuild_drop(0,rgx) 
      
! Define fix Z coordinate shift     

      
      fix_shift = boundary(3)/2.

      print '(/a,f10.7/)', "  * Constant shift in drop profile in z of ",fix_shift

#endif /*system=1: droplet*/

#endif /*profiles*/

#if STORE == 1 /* I/O refold  */
    print '(/a/)', "  * Storing unfolded coordinates "
!
!  Read unfold conf
!
        inquire(file='conf_unfold',exist=ufile )
        if(ufile) then    ! if conf_unfold exists, read it  
        print '(/a/)' ,'  *  Reading file conf_unfold '
        open (unit=180,file='conf_unfold')
        do i_part = 1 , n_part
        read(180,'(3e16.8)') r0_unfold(:,i_part)
        end do
        close (180)
        else ! if conf_unfold does not exist, generate the first  unfolded configuration 

        print '(/a/)' ,'  *  conf_unfold not present in dir. Unfolding first conf and starting from here '

! Init refold 
    

#if SYMMETRY == 0
        call refold(0,n_mon,1,part_init_d) 
        call refold(1,n_mon,1,part_init_d) 
        call refold(1,n_mon_d,1+part_init_d,part_init_e)
#elif SYMMETRY == 1 /*bulk */
        call refold(0,n_mon_d,1,part_init_e) 
        call refold(1,n_mon_d,1+part_init_d,part_init_e)
#endif

#ifdef PARTICLE_4
        call refold(1,n_mon_e,1+part_init_e,n_mon_tot) 
#endif
        
        endif
#endif

#if SYMMETRY == 1 

! Pressure tensor set to zero. Added in each routine that computes forces.
! See viscosity.f90 

      press_tensor(:,:) = 0. 
      press_tensor_mean(:,:) = 0. 
#endif 

#if SYSTEM == 0
#       ifdef PARTICLE_4
!        Mean velocity of particle 4 

               v4_mean(:) = 0.0 
               v4_mean2(:) = 0.0 
                 c4 = 0 

#       endif
#endif
#        ifdef DIFF
!       
!        Init diffusion Routine (only for particle 4 )
!       

             call diff_coef(1,r_time)
#       endif             
 
      end subroutine init_obser
