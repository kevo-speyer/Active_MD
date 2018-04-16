!    Here some physical quantities are calculated  
      subroutine observation
#include  'control_simulation.h'
      use commons
      use functions ! CHANGED BY KEVO TO RUN DROPPLET MODE 
! OBSOLETE #ifndef  PROFILES    
! OBSOLETE       use util
! OBSOLETE #else      
#ifdef PROFILES 
      use functions
#endif      
      implicit none
      integer :: i_chain_m,i
      real(kind=8) :: inv_count_obs
#if SYSTEM==1
  real(kind=8) :: v_cm_curr(3),r_trasl(3),curr_r_cm(3),rgx
#endif
!
      count_obs = count_obs + 1
!
      time_ave_count = time_ave_count + 1
!
!      if(time_ave_count.ge.n_time_ave) then
!
!!obs       if(f_logscale.eq.1) then
!!obs        n_time_ave=1.1*i_time
!!obs        if(n_time_ave.gt.200) then
!!obs         n_time_ave=200
!!obs         f_logscale=0
!!obs        end if
!!obs       end if
!!obs       if(f_logscale.eq.0) then
!!obs        time_ave_count = 0
!!obs       end if
!
!     periodic boundary refolding (to measure correctly the distances) )
!
!WARN: this should be translated to multiple chain length !!! 

! Brushes
#if SYMMETRY ==0 
        do i_part=1,part_init_d
            do i_dim=1,n_dim
                spabs(i_dim,i_part)=r0(i_dim,i_part)
                if(mod(i_part-1,n_mon).ne.0) then
                    2050      continue
                    if(spabs(i_dim,i_part)-spabs(i_dim,i_part-1).gt.               &
                        &half_boundary(i_dim)) then
                    spabs(i_dim,i_part)=spabs(i_dim,i_part)-boundary(i_dim)
                    goto 2050
                end if
                if(spabs(i_dim,i_part-1)-spabs(i_dim,i_part).gt.              &
                    &half_boundary(i_dim)) then
                spabs(i_dim,i_part)=spabs(i_dim,i_part)+boundary(i_dim)
                goto 2050
            end if
        end if
    end do
end do
#endif 

!
!      Angle of tilting and Re^2
!
!! obsoleted by now       if (f_angle.eq.1) then
!! obsoleted by now        angle=0.0
!! obsoleted by now        cntr=0
!! obsoleted by now        rendtot=0.0
!! obsoleted by now        zcoord=0.0
!! obsoleted by now        zcoord2=0.0
!! obsoleted by now       do i_part=1,n_mon*n_chain/2
!! obsoleted by now        if (mod(i_part-1,n_mon).eq.0) then
!! obsoleted by now         cntr=cntr+1
!! obsoleted by now         do i_dim=1,3
!! obsoleted by now          crdna(i_dim)=spabs(i_dim,i_part)
!! obsoleted by now          crdnb(i_dim)=spabs(i_dim,i_part+n_mon-1)
!! obsoleted by now         end do
!! obsoleted by now          rend(1)=crdnb(1)-crdna(1)
!! obsoleted by now          rend(2)=crdnb(2)-crdna(2)
!! obsoleted by now          rend(3)=crdnb(3)-crdna(3)
!! obsoleted by now          rendl2=rend(1)**2+rend(2)**2+rend(3)**2
!! obsoleted by now          zcomp=crdna(3)-crdnb(3)
!! obsoleted by now          xcomp=crdna(1)-crdnb(1)
!! obsoleted by now!
!! obsoleted by now!
!! obsoleted by now        if(xcomp.lt.0) then
!! obsoleted by now          angle=angle-180*atan(zcomp/xcomp)/3.141592
!! obsoleted by now        end if
!! obsoleted by now        if(xcomp.gt.0) then
!! obsoleted by now          angle=angle+180-180*atan(zcomp/xcomp)/3.141592
!! obsoleted by now        end if
!! obsoleted by now        if(xcomp.eq.0) then
!! obsoleted by now          angle=angle+90
!! obsoleted by now        end if
!! obsoleted by now          rendtot=rendtot+rendl2
!! obsoleted by now          zcoord=zcoord+rend(3)
!! obsoleted by now        end if
!! obsoleted by now       end do
!! obsoleted by now         do i_part=n_mon*n_chain/2+1,n_mon*n_chain
!! obsoleted by now         if (mod(i_part-1,n_mon).eq.0) then
!! obsoleted by now         cntr2=cntr2+1
!! obsoleted by now         do i_dim=1,3
!! obsoleted by now          crdna2(i_dim)=spabs(i_dim,i_part)
!! obsoleted by now          crdnb2(i_dim)=spabs(i_dim,i_part+n_mon-1)
!! obsoleted by now         end do
!! obsoleted by now          rend2(1)=crdnb2(1)-crdna2(1)
!! obsoleted by now          rend2(2)=crdnb2(2)-crdna2(2)
!! obsoleted by now          rend2(3)=crdnb2(3)-crdna2(3)
!! obsoleted by now          rendl22=rend2(1)**2+rend2(2)**2+rend2(3)**2
!! obsoleted by now          zcomp2=crdnb2(3)-crdna2(3)
!! obsoleted by now          xcomp2=crdna2(1)-crdnb2(1)
!! obsoleted by now!
!! obsoleted by now        if(xcomp2.lt.0) then
!! obsoleted by now          angle2=angle2-180*atan(zcomp2/xcomp2)/3.141592
!! obsoleted by now        end if
!! obsoleted by now        if(xcomp2.gt.0) then
!! obsoleted by now          angle2=angle2+180-180*atan(zcomp2/xcomp2)/3.141592
!! obsoleted by now        end if
!! obsoleted by now        if(xcomp2.eq.0) then
!! obsoleted by now          angle2=angle2+90
!! obsoleted by now        end if
!! obsoleted by now          rendtot2=rendtot2+rendl22
!! obsoleted by now          zcoord2=zcoord2+rend2(3)
!! obsoleted by now        end if
!! obsoleted by now       end do
!! obsoleted by now       end if



!     force acting on top wall with explicit wall atoms 
!!! obsolete by now  if(f_explicit_wall) then !if not explicit wall other values are calculated behind
!!! obsolete by now        if(f_force.eq.1) then
!!! obsolete by now          write(80,*) r_dummy,ftw(1),ftw(2),ftw(3)
!!! obsolete by now        end if
!!! obsolete by now  end if       
!!! obsolete by now      end if ! time>time_ave
!
!
!
      v_wall_wall_1 = v_wall_wall_1 + v_wall_wall
      t_wall_1      = t_wall_1      + t_wall
      e_wall_wall_1 = e_wall_wall_1 + (v_wall_wall+t_wall)
!
      v_wall_wall_2 = v_wall_wall_2 + v_wall_wall**2
      t_wall_2      = t_wall_2      + t_wall**2
      e_wall_wall_2 = e_wall_wall_2 + (v_wall_wall+t_wall)**2
!
      v_fluid_fluid_1 = v_fluid_fluid_1 +                               &
     &(v_fluid_fluid + v_intra_molec + v_fluid_wall)
      t_fluid_1      = t_fluid_1      + t_fluid
      e_fluid_fluid_1 = e_fluid_fluid_1 + (v_fluid_fluid+t_fluid)
!
      v_fluid_fluid_2 = v_fluid_fluid_2 +                               &
     &(v_fluid_fluid + v_intra_molec + v_fluid_wall)**2
      t_fluid_2      = t_fluid_2      + t_fluid**2
      e_fluid_fluid_2 = e_fluid_fluid_2 + (v_fluid_fluid+t_fluid)**2
!
      v_total_1 = v_total_1 + v_total
      t_total_1 = t_total_1 + t_total
      e_total_1 = e_total_1 + e_total
!
      v_total_2 = v_total_2 + v_total**2
      t_total_2 = t_total_2 + t_total**2
      e_total_2 = e_total_2 + e_total**2
!
      z_1 = z_1 + r0_twall(n_dim)
      z_2 = z_2 + r0_twall(n_dim)**2
!claudio
 !      print *,v_intra_molec,v_intra_molec_d
      v_tot_intra_1 = v_tot_intra_1 + v_intra_molec
      v_tot_intra_2 = v_tot_intra_2 + v_intra_molec**2
      v_tot_intra_d_1 = v_tot_intra_d_1 + v_intra_molec_d
      v_tot_intra_d_2 = v_tot_intra_d_2 + v_intra_molec_d**2
!      print*,v_tot_intra_1,v_tot_intra_d_1
!fin claudio

! -------------------------------------------------------------------------------      

!TO FIX: to multiple chain lengths
! This does not calculates correctly for  different chain types and  lengths

    spabs(:,:) = r0(:,:)

#if SYMMETRY == 0       /* if channel geom and brushes*/ 
            do i_part=1,part_init_d
                do i_dim=1,n_dim
                    if(mod(i_part-1,n_mon).ne.0) then
                        if(spabs(i_dim,i_part)-spabs(i_dim,i_part-1).gt.half_boundary(i_dim)) then
                            spabs(i_dim,i_part)=spabs(i_dim,i_part)-boundary(i_dim)
                            cycle
                        end if
                        if(spabs(i_dim,i_part-1)-spabs(i_dim,i_part).gt.half_boundary(i_dim)) then
                            spabs(i_dim,i_part)=spabs(i_dim,i_part)+boundary(i_dim)
                            cycle
                        end if
                    end if
                end do
            end do
#endif     
        do i_part=part_init_d+1,part_init_e
            do i_dim=1,n_dim
                if(mod(i_part-1,n_mon_d).ne.0) then
                    if(spabs(i_dim,i_part)-spabs(i_dim,i_part-1).gt.half_boundary(i_dim)) then
                        spabs(i_dim,i_part)=spabs(i_dim,i_part)-boundary(i_dim)
                        cycle
                    end if
                    if(spabs(i_dim,i_part-1)-spabs(i_dim,i_part).gt.half_boundary(i_dim)) then
                        spabs(i_dim,i_part)=spabs(i_dim,i_part)+boundary(i_dim)
                        cycle
                    end if
                end if
            end do
        end do
#ifdef PARTICLE_4
        do i_part=part_init_e+1,n_mon_tot
            do i_dim=1,n_dim
                if(mod(i_part-1,n_mon_e).ne.0) then

                    if( spabs(i_dim,i_part)-spabs(i_dim,i_part-1) >  half_boundary(i_dim) ) then
                        spabs(i_dim,i_part)=spabs(i_dim,i_part)-boundary(i_dim)
                        cycle
                    end if

                    if(spabs(i_dim,i_part-1)-spabs(i_dim,i_part).gt.half_boundary(i_dim) ) then
                        spabs(i_dim,i_part)=spabs(i_dim,i_part)+boundary(i_dim)
                        cycle
                    end if

                end if
            end do
        end do
#endif

    !end do ! ndim

!-- end refolding

#if SYSTEM == 0 || SYSTEM == 3
!
! Done only for brushes 
!

    r_g2(:,:) = 0.d0
!
! **** Rg for Brushes
!
 
        do i_chain = 1,n_chain
            v_dummy(:) = 0.0
            r_cm(:) = sum(spabs(:,(i_chain-1)*n_mon+1:i_chain*n_mon),dim=2)*inv_n_mon !  / dble(n_mon)

            do i_part=n_mon*(i_chain-1)+1 ,i_chain*n_mon
                v_dummy(:) = v_dummy(:) +   (spabs(:,i_part) - r_cm(:))**2 
            end do
            v_dummy(:) = v_dummy(:)*inv_n_mon  ! Rg^2 of one polymer
            r_g2(1,:) =  r_g2(1,:) + v_dummy(:)  ! mean Rg*N
        end do
        r_g2(1,:) =  r_g2(1,:)* inv_n_chain ! /dble(n_chain)! mean Rg
!
! **** Rg for Melt
!
        do i_chain = n_chain+1,n_chain + n_chain_d
            v_dummy(:) = 0.0
            r_cm(:) = sum(spabs(:,(i_chain-1)*n_mon_d+1:i_chain*n_mon_d),dim=2) / dble(n_mon_d)
            do i_part= (i_chain-1)*n_mon_d+1,i_chain*n_mon_d
                !     v_dummy(:) = v_dummy(:) +   spabs(:,i_part)**2 - r_cm(:)**2   ! WRONG 
                v_dummy(:) = v_dummy(:) +   (spabs(:,i_part) - r_cm(:))**2  
            end do
            v_dummy(:) = v_dummy(:) / dble(n_mon_d)! Rg of one polymer
            r_g2(2,:) =  r_g2(2,:) + v_dummy(:)! sumN Rg
        end do
        r_g2(2,:) =  r_g2(2,:)/dble(n_chain_d) ! mean Rg^2

        r_g2_mean(1,:) = r_g2_mean(1,:)+ r_g2(1,:) !brushes
        r_g2_mean(2,:) = r_g2_mean(2,:)+ r_g2(2,:) ! droplet/melt
        r_g2_mean(3,:) = r_g2_mean(3,:)+ (r_g2(1,:) + r_g2(2,:))/2.0 ! total 

        r_par_per2(1) = r_par_per2(1)+ (r_g2(1,3)- 0.5*(r_g2(1,1)+r_g2(1,2)))/sum(r_g2(1,:))
        r_par_per2(2) = r_par_per2(2)+ (r_g2(2,3)- 0.5*(r_g2(2,1)+r_g2(2,2)))/sum(r_g2(2,:))

#endif /* SYSTEM=0 or system=3 */
! ------------ end gyration stuff

! energies out

! Total momentum of the system
        vec_dummy(:) = 0.
        do i_part = 1 , n_mon_tot 
            vec_dummy(:) = vec_dummy(:) + mass(i_part)*v(:,i_part)
        end do
        tot_P(:) = tot_P(:) + vec_dummy(:)

! Drop data to disk 

        if(mod(i_time-n_relax,5).eq.0) then

            inv_count_obs =1./dble(count_obs)

            write(40,'(i7,6f15.5)')i_time,                       &
                v_total_1 *inv_count_obs , & 
                t_total_1*inv_count_obs  , &  
                e_total_1*inv_count_obs  , &
                tot_P(:)*inv_count_obs
            write(41,'(i7,3f15.5)')i_time,           &
                v_fluid_fluid ,  & 
                v_intra_molec ,  & 
                v_fluid_wall  

            write(42,'(i7,8f15.7)')i_time,                       &
                r_g2_mean(1,:)  *inv_count_obs           ,      & !gyration vector Rg brushes
                sum(r_g2_mean(1,:))  *inv_count_obs           , & ! brush modulus of Rg^2 
                r_g2_mean(2,:) *inv_count_obs            ,      &  !gyr vector drop/melt
                sum(r_g2_mean(2,:))  *inv_count_obs               ! melt/droplet modulus of Rg^2 
            write(43,'(i7,2f15.7)')i_time,r_par_per2(1:2) *inv_count_obs 

!---  force on walls if not explicit wall atoms

#       if WALL == 2 && WR_FORCES == 0 /* implicit wall LJ 9-3 and  Writes sum of forces on brush heads*/

!    if(.not.f_explicit_wall) then
    
!  ---  Write  force on _heads. They are taken in  corrector.f90 or verlet_velocities.f90

write(80,'(7e18.7)') r_time,                                             &
                     -sum(f_on_heads(:,1:n_chain_2) ,dim=2) ,            & ! upper wall
                     -sum(f_on_heads(:,n_chain_2+1:n_chain),dim=2)         ! lower wall


! Modified by Kevo 4/2014 write force on each head in fort.80 
! Total force on the brush heads (negated). Separate in lower and upper walls
#       endif
#       if WALL == 2 && WR_FORCES == 1 /* implicit wall LJ 9-3 and Writes individual force on brush heads*/

    write(80,*) n_chain_2
    do i_part=1,n_chain_2
        write(80,'(6e18.7)') f_on_heads(:,i_part)   ,           &  ! upper wall
                         f_on_heads(:,n_chain_2+i_part)         ! lower wall
    end do                   
                         !!    end if !not explicit wall
#       endif

end if ! mod(i_time-n_relax,5).eq.0

! ---------- Profiles Calculation 

#if SYSTEM == 0 || SYSTEM == 2 || SYSTEM == 3 /*channel or charges*/
#       ifdef PROFILES 
! Densities (order N algorithm )


        call histo(2,dim_prof_dens) 

! Velocities  (order N algorithm )


       call  velocity_prof(2,dim_prof_dens,zz_min,boundary(3)) 

#       endif
!
! Measurements of quantities of particle 4 
!
#   ifdef PARTICLE_4
#       if SYSTEM == 0
          
             do i_part=1,n_tot_e            
                   c4 = c4 + 1 

#           ifdef PROFILES
! if v is defined 
                    v4_mean(:) = v4_mean(:) + v(:,part_init_e+1+i_part) 
                    v4_mean2(:) = v4_mean2(:) + v(:,part_init_e+1+i_part)**2
#           else
                    v4_mean(:) = v4_mean(:) + v(:,part_init_e+1+i_part) 
                    v4_mean2(:) = v4_mean2(:) + v(:,part_init_e+1+i_part)**2
#           endif /* profiles */
             end do  

#       endif
#    endif /* particle 4 */

#       ifdef DIFF
!        Diffusion of particle 4 
                   call diff_coef(3,r_time)
#       endif /* DIFF */            

#endif ! system=0: channel

#if SYSTEM == 1 /* droplets */

! --- DROP CM Velocity and position  calculation

! NOTE: assumes that the all the free chains belong to the droplet 
! (erroneous if there is significant vapor phase)

           v_cm_curr(:) = 0.
           do i_part = part_init_d+1,n_mon_tot !loop over melt particles
               v_cm_curr(:) = v_cm_curr(:) + v(:,i_part)
           end do
           
           v_cm_curr(:) = v_cm_curr(:) / dble(nm) 
                !ori:                  / dble(part_init_d), bad normalization !? 
           
           vcm_d(:) = vcm_d(:) + v_cm_curr(:)
           vcm_d2(:) = vcm_d2(:) + v_cm_curr(:)**2



! Write out droplet's CM position and velocity 
! rebuild_drop2(2)  should be used. It's fast and gives a continuous trajectory
! for droplet
!        call rebuild_drop2(1) ! original. drop_cm(:) and n_cm are recalculated !obsolete
                                    !Calculations done beeter in rebuild_drop2(2)

        call rebuild_drop2(2) ! More accuratelly and efficient calculation of drop_cm(:) and n_cm
                        
!        call rebuild_drop2(3) ! Another way to calculate of drop_cm(:) and n_cm
!        call rebuild_drop2(4) ! Another way to calculate of drop_cm(:) and  n_cm



    
        write(54,'(i7,6(g16.8,2x))') i_time,drop_cm(:),v_cm_curr(:)

! Droplet profiles [ Taken from mide.f90 for mide_drop_force]

#ifdef PROFILES 
! ------- Unfold chains ( each chain is rebuilded removing if necessary the PBC )
! And defined for the first time r0_unfold

!ori mfa_drop_force        call refold(2) 

        call refold(1,n_mon_d,1+part_init_d,n_mon_tot) !get r0_unfold right for melt

        call refold(1,n_mon,1,part_init_d) !Add by Kevo 5/2016, get r0_unfold right for brush


! ------- Calculates all the CMs of the chains 

        call calc_cms(2) !DEFINE LEGACY


! ------- Assemble the chains according to the relative distance of the drop CM.
!         r0_unfold recalculated.   Also calculates Rg_drop()

        call  rebuild_drop(n_cm,rgx)

! Calculates and write out Rg_z and Rg_x

        call calc_drop_rg(rgx)

! The substraction of CM is done only in X 

!only in X        r_trasl(1) = drop_cm(1) !old rcm(r0_unfold(1:3,n_mon*n_chain+1:n_part))
!only in X        r_trasl(2:3) = 0.

! The substraction of CM is done only in X AND Z
! NOTE: a fix shift in z is done in order to force positive Z coordinate for all the particles. 

        r_trasl(1) = drop_cm(1) - boundary(1)/2. !Added: - boundary(1)/2.
        r_trasl(2) = 0.
        r_trasl(3) = drop_cm(3) - fix_shift  !fix_shift = boundary(3)/2. OK 

!
! --- Translation of coordinates so that drop  CM  in (/ Lx/2,0,Lz/2 /)
!
!OBSOLETE?        r_trasl(:)  = r_trasl(:) - Lx_o_2(1:3) 
!OBSOLETE?
!OBSOLETE?        r0_unfold(:,:) = r0_unfold(:,:) - spread ( r_trasl(1:3),1,n_part ) ! all

!WARNING RESET r0_unfold to get a good histogram
    r0_unfold(:,:) = r0(:,:)

! Melt translation
! Translates in X and Z

        do i=part_init_d+1,n_mon_tot
            r0_unfold(:,i) = r0_unfold(:,i) - r_trasl(:) 
#ifndef NO_WARNS            
            if (r0_unfold(3,i) <0.) then
                print *,i,r0_unfold(3,i),"PROFS: out of BOX in Z"
            end if
            if (r0_unfold(3,i) >boundary(3) ) then
                print *,i,r0_unfold(3,i),"PROFS: OUT of BOX in Z. Taking it back "
                r0_unfold(:,i) = boundary(3) -1.2
            end if
#endif            
        end do

! Brush translation
! Translated  only in X. 
! No translation done for brush, because it gives a blurry image

!        do i=1,part_init_d
!             r0_unfold(1,i) = r0_unfold(1,i) - r_trasl(1)
!        end do

! Translate X positions the middle of the box. Needed after drop rebuilding.

!     do i = 1 , n_part
!        r0_unfold(:,i) = r0_unfold(:,i) !+  lx_o_2(1:3) 
!     end do
!
! Density prof calculation
!
!  Refold in y before getting the density 

        call  coor_refold(2) !Apply PBC in y 
!        call  coor_refold(1) !Apply PBC for a box length = 2 * boundary(1)
        call  coor_refold(4) !Apply PBC in x for a box length =  boundary(1)

! ------ Density profile 
! Order N algorithm: 

         call dens_prof(3,dens_xz(:,:,:),dens_xz_b(:,:,:),dens_xz_2(:,:,:),dens_xz_b_2(:,:,:),n_box_xz,r_box_xz) 

! ----- Velocity profile


          call vel_prof_2d(3,vel_xz,vel_xz2,vel_xz_b,vel_xz2,n_box_xz,r_box_xz,dens_xz) ! 

#       endif /* profiles */
#endif /* droplets */
#       if SYMMETRY == 1 
          ! Calculates viscosity, pressure and presure tensor
          call viscosity()
#endif

 end subroutine observation
