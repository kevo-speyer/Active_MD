      subroutine obser_out
#include 'control_simulation.h'     
      use commons 
      use ziggurat
! OBSOLETE #ifndef PROFILES      
! OBSOLETE       use util , only: velocity_prof
! OBSOLETE #else      
#ifdef PROFILES
      use functions
#endif      
      implicit none
      real(kind=8) :: viscosity ,inv_count_obs 
#ifdef PROFILES
     real(kind=8) :: vol_layer,dz 
     real (kind=8) , allocatable :: sig2_v_prof(:,:,:)
     integer :: i
#endif
     real(kind=8) :: v4m(3),v4_2m(3)


      write(20,*)
      write(20,220) n_fcell_max,n_bin_fl,                               &
       "  max. observed and maximum allowed fl. binnings"
       write(20,220) n_wcell_max,n_bin_wa,                               &
       "  max. observed and maximum allowed wa. binnings"
       write(20,*)
       write(20,220) n_neigh_fl_max,n_neigh_fl,                          &
       "  max. observed and maximum allowed fl fl neighbors"
       write(20,220) n_neigh_fw_max,n_neigh_wa,                          &
       "  max. observed and maximum allowed fl wa neighbors"
       write(20,220) n_neigh_ww_max,n_neigh_ww,                          &
       "  max. observed and maximum allowed wa wa neighbors"
!
!if(f_explicit_wall) then   !only if considering explicit wholw atoms
#   if WALL == 1 /* explicit wall */
      v_wall_wall_1 = v_wall_wall_1/n_wall/count_obs
      t_wall_1      = t_wall_1     /n_wall/count_obs
      e_wall_wall_1 = e_wall_wall_1/n_wall/count_obs
!
      v_wall_wall_2 = v_wall_wall_2/n_wall**2/count_obs
      t_wall_2      = t_wall_2     /n_wall**2/count_obs
      e_wall_wall_2 = e_wall_wall_2/n_wall**2/count_obs
!
      v_wall_wall_2 = v_wall_wall_2 - v_wall_wall_1**2
      t_wall_2      = t_wall_2      - t_wall_1**2
      e_wall_wall_2 = e_wall_wall_2 - e_wall_wall_1**2
!
      v_wall_wall_2 = n_wall*v_wall_wall_2/temp**2
      t_wall_2      = n_wall*t_wall_2     /temp**2
      e_wall_wall_2 = n_wall*e_wall_wall_2/temp**2
#  else ! if no explicit wall atoms. We take in a different form the mean values
!warning: this should be check
 
inv_count_obs = 1./count_obs

      v_wall_wall_1 = v_wall_wall_1*inv_count_obs
      t_wall_1      = t_wall_1     *inv_count_obs
      e_wall_wall_1 = e_wall_wall_1*inv_count_obs
!
      v_wall_wall_2 = v_wall_wall_2**2*inv_count_obs
      t_wall_2      = t_wall_2     **2*inv_count_obs
      e_wall_wall_2 = e_wall_wall_2**2*inv_count_obs
!
      v_wall_wall_2 = v_wall_wall_2 - v_wall_wall_1**2
      t_wall_2      = t_wall_2      - t_wall_1**2
      e_wall_wall_2 = e_wall_wall_2 - e_wall_wall_1**2
!
      v_wall_wall_2 = v_wall_wall_2/temp**2
      t_wall_2      = t_wall_2     /temp**2
      e_wall_wall_2 = e_wall_wall_2/temp**2
!
#   endif ! explicit wall atoms

      write(20,*)
      write(20,202) v_wall_wall_1,v_wall_wall_2," wall_wall potential"
      write(20,202) t_wall_1     ,t_wall_2     ," wall_wall kinetic e"
      write(20,202) e_wall_wall_1,e_wall_wall_2," wall_wall total ene"
!
      v_fluid_fluid_1 = v_fluid_fluid_1/n_mon_tot/count_obs
      t_fluid_1       = t_fluid_1      /n_mon_tot/count_obs
      e_fluid_fluid_1 = e_fluid_fluid_1/n_mon_tot/count_obs
!
      v_fluid_fluid_2 = v_fluid_fluid_2/n_mon_tot**2/count_obs
      t_fluid_2       = t_fluid_2      /n_mon_tot**2/count_obs
      e_fluid_fluid_2 = e_fluid_fluid_2/n_mon_tot**2/count_obs
!
      v_fluid_fluid_2 = v_fluid_fluid_2 - v_fluid_fluid_1**2
      t_fluid_2       = t_fluid_2       - t_fluid_1**2
      e_fluid_fluid_2 = e_fluid_fluid_2 - e_fluid_fluid_1**2
!
      v_fluid_fluid_2 = n_mon_tot*v_fluid_fluid_2/temp**2
      t_fluid_2       = n_mon_tot*t_fluid_2      /temp**2
      e_fluid_fluid_2 = n_mon_tot*e_fluid_fluid_2/temp**2
!clau intramolecular
      v_tot_intra_1 = v_tot_intra_1 /n_mon_tot/count_obs

    if(n_mon_d>0.and.n_chain_d>0) then
          v_tot_intra_d_1 = v_tot_intra_d_1/(dble(n_mon_d*n_chain_d))/count_obs
          v_tot_intra_d_2 = v_tot_intra_d_2 /(n_mon_d*n_chain_d)**2/count_obs    
    else ! if there are no monomers, v_intras to 0 
          v_tot_intra_d_1 = 0.  
          v_tot_intra_d_2 = 0. 
    end if
      v_tot_intra_2 = v_tot_intra_2     /n_mon_tot**2/count_obs                  
      v_tot_intra_2 = v_tot_intra_2 - v_tot_intra_1**2   
      v_tot_intra_d_2 = v_tot_intra_d_2 - v_tot_intra_d_1**2
      v_tot_intra_2 = v_tot_intra_2      /temp**2
      v_tot_intra_d_2 = v_tot_intra_d_2  /temp**2
! warn: as it is calculated per particle. Both contributions
! should not be very different
      write(20,*)
      write(20,'(2f13.4,a40)') v_tot_intra_1,v_tot_intra_2, "total intramolecular energy per part."
      write(20,'(2f13.4,a50)') v_tot_intra_d_1,v_tot_intra_d_2, "intramolecular energy droplet/melt per part"
      write(20,*)
      write(20,202) v_fluid_fluid_1,v_fluid_fluid_2," fl_fl potential"
      write(20,202) t_fluid_1      ,t_fluid_2      ," fl_fl kinetic e"
      write(20,202) e_fluid_fluid_1,e_fluid_fluid_2," fl_fl total ene"
!
      v_total_1 = v_total_1/n_part/count_obs
      t_total_1 = t_total_1/n_part/count_obs
      e_total_1 = e_total_1/n_part/count_obs
!
      v_total_2 = v_total_2/n_part**2/count_obs
      t_total_2 = t_total_2/n_part**2/count_obs
      e_total_2 = e_total_2/n_part**2/count_obs
!
      v_total_2 = v_total_2 - v_total_1**2
      t_total_2 = t_total_2 - t_total_1**2
      e_total_2 = e_total_2 - e_total_1**2
!
      v_total_2 = n_part*v_total_2/temp**2
      t_total_2 = n_part*t_total_2/temp**2
      e_total_2 = n_part*e_total_2/temp**2
!
      write(20,*)
      write(20,202) v_total_1,v_total_2," total potential"
      write(20,202) t_total_1,t_total_2," total kinetic e"
      write(20,202) e_total_1,e_total_2," total total ene"
!
      z_1 = z_1/count_obs
      z_2 = z_2/count_obs
      z_2 = z_2 - z_1**2
      z_2 = (z_2/temp)*n_part
      write(20,*)
      write(20,202) z_1,z_2," average spacing"
!
      write(20,*)
!!!ori      write(20,260) i_obs_l1,i_obs_l2,i_obs_l3,i_obs_l4,i_obs_l5,       &
!!!ori     &              i_obs_l6,"  diffusive observations"
! Radious of gyration stuff
  write(20,*) !space 
  write(20,'(3f15.5,a60)') r_g2_mean(1,:)  /dble(count_obs) , "Mean Rg^2 for brushes"          
  write(20,'(3f15.5,a60)')r_g2_mean(2,:) /dble(count_obs)   , "Mean Rg^2 for droplet/melt"          
  write(20,'(3f15.5,a60)')r_g2_mean(3,:) /dble(count_obs)   , "Mean total Rg^2 "
  write(20,'(f15.5,a60)') r_par_per2(1)/dble(count_obs), "Mean Rg(z)^2-0,5*(Rg^2(x)+Rg^2(y))/sum(Rg^2) brushes"          
  write(20,'(f15.5,a60)') r_par_per2(2)/dble(count_obs), "Mean Rg(z)^2-0,5*(Rg^2(x)+Rg^2(y))/sum(Rg^2) melt"          
! Force on heads if there is no explicit wall

  write(20,*) !space 
!if(.not.f_explicit_wall) then
#       if WALL != 1 
         f_on_heads(:,:) = f_on_heads(:,:)/dble(count_obs)
         if(n_chain>0) then ! If there are brushes
             mean_f_on_heads(:) = abs(sum(f_on_heads(:,1:n_chain/2))) & ! top wall
                 + abs( sum(f_on_heads(:,n_chain/2+1:n_chain),dim =1) )/2. ! botom wall

             write(20,'(3f16.5,a40)') mean_f_on_heads(:), "Mean total force on the wall  "
! ---   shear stress assuming v=vx     
             write(20,'(f17.7,a40)') (mean_f_on_heads(1)+mean_f_on_heads(3))/(2.*surface), " Shear stress "
         end if
         !    this is true assuming the velocity is in x direction
         write(20,*) !space 
         if(va_spring_twall(1) /= 0.) then
             viscosity = z_space_wall* (mean_f_on_heads(1))/(2.*surface*va_spring_twall(1)) 
             write(20,'(f16.5,a40)')  viscosity," Effective viscocity "
         end if
#        endif
  write(20,*) !space 
#if SYSTEM == 0 
#      ifdef PARTICLE_4

! Measurements of particle 4

       v4m = v4_mean(:)/dble(c4)
       v4_2m = v4_mean2(:)/dble(c4)    
      write(20,'(/8x,a/)') "Quantities for particle 4" 
      write(20,'(8x,a,3g15.6,i8)') , "<v4>,N =",v4m,c4
      write(20,'(8x,a,3g15.6)') ,    "<v4^2> =",v4_2m
      write(20,'(8x,a,3g15.6)') ,    "sig(v4)=",sqrt(v4_2m(:)-v4m(:)**2)
      write(20,'(8x,a,3g15.6)') ,    "err(v4)=",sqrt(v4_2m(:)-v4m(:)**2)/sqrt(dble(c4))
#       endif /* PARTICLE 4 */
#endif
! *** Write Velocity Profile 

#ifdef PROFILES
!
! Density profile 
!
#if SYSTEM == 0 || SYSTEM == 2  || SYSTEM == 3 /*channel or charges */

       call  histo(1,dim_prof_dens)
!
! Vel profile
!
! NOTE: this call MUST be done after histo(1 ,...)

       call  velocity_prof(1,dim_prof_dens,zz_min,boundary(3) ) !z_space_wall,v_prof(:,:,:),v_prof_2(:,:,:))

! 
!  ---- Normalized velocity profile  and  sig2 vel profile   
! 
#ifdef PARTICLE_4
       allocate (sig2_v_prof(dim_prof_dens,3,3))
#else
       allocate (sig2_v_prof(dim_prof_dens,3,2))
#endif

       vol_layer = surface*boundary(3)/dble(dim_prof_dens)

       open (unit=75,file='norm_vel_prof.mide',status='unknown')
       open (unit=77,file='norm_sig2_v_prof.mide',status='unknown')


! NOTE: in order to get a "temperature profile we have : T=mass*sig2_v. 
! We multiply by  the mass 

       do i = 1,dim_prof_dens
!       
! --- Normalization of sigma and velocity profiles
!
           if (histo_out(i,1) >0.) then
              v_prof(i,:,1) = v_prof(i,:,1)/histo_out(i,1)/ vol_layer  ! brush 
              v_prof_2(i,:,1) = v_prof_2(i,:,1)/histo_out(i,1)/ vol_layer  ! brush 
              sig2_v_prof(i,:,1) = (v_prof_2(i,:,1) - v_prof(i,:,1)**2) 
           else  
              v_prof(i,:,1) = 0.
              sig2_v_prof(i,:,1) = 0. 
           end if

       if (histo_out(i,2) >0.) then
          v_prof(i,:,2) = v_prof(i,:,2)/histo_out(i,2) / vol_layer ! melt
          v_prof_2(i,:,2) = v_prof_2(i,:,2)/histo_out(i,2) / vol_layer ! melt
          sig2_v_prof(i,:,2) = (v_prof_2(i,:,2) - v_prof(i,:,2)**2) 
      else  
          v_prof(i,:,2) = 0.
          sig2_v_prof(i,:,2) = 0. 
      end if

#ifdef PARTICLE_4
       if (histo_out(i,3) >0.) then
          v_prof(i,:,3) = v_prof(i,:,3)/histo_out(i,3) / vol_layer ! part4
          v_prof_2(i,:,3) = v_prof_2(i,:,3)/histo_out(i,3) / vol_layer ! part4
          sig2_v_prof(i,:,3) = (v_prof_2(i,:,3) - v_prof(i,:,3)**2) 
      else  
          v_prof(i,:,3) = 0.
          sig2_v_prof(i,:,3) = 0. 
      end if
#endif

       dz = zz_min +  boundary(3)*(dble(i-1)+0.5)/dble(dim_prof_dens)

! Normalize with the masses . The output will be T(z)

       sig2_v_prof(i,:,1) = sig2_v_prof(i,:,1)*mass_type(2) ! brush
       sig2_v_prof(i,:,2) = sig2_v_prof(i,:,2)*mass_type(3) ! melt
#ifdef PARTICLE_4       
       sig2_v_prof(i,:,3) = sig2_v_prof(i,:,3)*mass_type(4) ! particle 4
#endif

#ifndef PARTICLE_4
! Write out vel prof
       write(75,'(7f16.10)') dz ,v_prof(i,1,1),v_prof(i,1,2),   &
                                 v_prof(i,2,1),v_prof(i,2,2),    &
                                 v_prof(i,3,1),v_prof(i,3,2)
! Write out temp prof
       write(77,'(7f16.10)') dz ,sig2_v_prof(i,1,1),sig2_v_prof(i,1,2),   &
                                 sig2_v_prof(i,2,1),sig2_v_prof(i,2,2),    &
                                 sig2_v_prof(i,3,1),sig2_v_prof(i,3,2)
#else
       write(75,'(10f16.10)') dz ,v_prof(i,1,1),v_prof(i,1,2), v_prof(i,1,3), &
                                  v_prof(i,2,1),v_prof(i,2,2), v_prof(i,2,3), &
                                  v_prof(i,3,1),v_prof(i,3,2), v_prof(i,3,3)

       write(77,'(10f16.10)') dz ,sig2_v_prof(i,1,1),sig2_v_prof(i,1,2),sig2_v_prof(i,1,3),   & ! X
                                  sig2_v_prof(i,2,1),sig2_v_prof(i,2,2),sig2_v_prof(i,2,3),    & ! Y
                                  sig2_v_prof(i,3,1),sig2_v_prof(i,3,2),sig2_v_prof(i,3,3)       ! Z
#endif
    end do ! over z bins
       close(75)
       close(77)
#endif /*SYSTEM=0: channel*/
#endif /*profiles*/

#   if SYSTEM==1
        
#       ifdef PROFILES 

     call dens_prof(1,dens_xz(:,:,:),dens_xz_b(:,:,:),dens_xz_2(:,:,:),dens_xz_b_2(:,:,:),n_box_xz, &
     r_box_xz)  

     call vel_prof_2d(1,vel_xz,vel_xz2,vel_xz_b,vel_xz2,n_box_xz,r_box_xz,dens_xz) ! 

#       endif /*profiles*/
#   endif /*system=1: droplets*/

# if SYMMETRY == 1 
! Pressure tensor stuff 
 
press_tensor_mean(:,:) = press_tensor_mean(:,:) / count_obs

write (20,'(/a/)') " ***** Mean pressure tensor components ******"

write(20,'(10x,3f16.7)') press_tensor_mean(1,:)
write(20,'(10x,3f16.7)') press_tensor_mean(2,:)
write(20,'(10x,3f16.7)') press_tensor_mean(3,:)

! Pressure given as tr(press_tensor) / 3.
write (20,'(/a,f16.7/)') " Mean pressure=", (press_tensor_mean(1,1) +press_tensor_mean(2,2) +press_tensor_mean(3,3) ) /3.
#endif 




! *** Closing storage  files 

      call store_config(5)

! Close a random number generator

        open(unit=222,file='random_seed.dat',status='unknown')
        write(222,*) shr3()

  202 format(2e15.6,a)
  220 format(2i8,a)

      end   subroutine obser_out
