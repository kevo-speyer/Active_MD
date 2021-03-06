subroutine calc_solv_solv_force()
#include 'control_simulation.h'
      use commons
#ifdef _OPENMP
      use omp_lib
      use Par_Zig_mod, only: par_rnor,par_uni
#else
      use ziggurat, only: rnor,uni
#endif

     implicit none
     real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc,f_ipart(3) ! Needed for DPD_EMBEDDED only
     real (kind=8) :: l_eps,r_cut2,r_61, sigma_2_ss
     real (kind=8) :: f_cou_real(3)
     real(kind=8)  :: inv_r_2,inv_sqrt_r_2
! cache blocking 
     integer :: ii_part,ii_neigh,i,j
      logical, parameter :: debug=.false.
!DEBUG
!print*,"Entro a calc_solv_solv_force"
!/DEBUG
!For solvent_solvent interactions only, so l_eps is fixed 
      l_eps = epsil(3,3)
      r_cut2 = range_2(3,3)
      sigma_2_ss = sigma_2(3,3)

     
#       if SYSTEM == 2 || SYSTEM == 3
            v_coul = 0.
#       endif

!!#ifdef DPD_EMBEDDED
! NOTE: this is at the beggining of dpd_forces_ll when the DPD forces are not
! embededd here. This must be added if dpd_forces is never called in the program
! flow. 

! Force switch-on update  

! 
!---  Define effective minimium distance used for force "switch on" at the beginning 
!     This is done for each time step

#ifdef FORCE_SWITCH_ON

      if(r_time.lt.r_2_min_time) then
       r_2_min = (1.-r_time/r_2_min_time)*r_2_min_init !*range_2(:,:)!ORI*r_2_min_init
      else
       r_2_min = 0.
      end if
#endif      
      v_sol_sol = 0.
!
! LJ fluid-fluid force and V  calculation 
!

!DEBUG

!BEGIN PARALLEL ZONE

!Warning: Paralelization not adapted for SYMMETRY=1 
# if BIN_TYPE == 0 || BIN_TYPE == 1
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(r_2_min,force_long,v_sol_sol,n_part,a_type,ff_list, range_2, r0,inv_boundary, boundary,sigma_2,e_shift,sig_long,mass,friction,v, inv_range_2, r_cut_ss, inv_r_cut_ss,r_cut_dpd_2)
#ifdef _OPENMP
    ith=omp_get_thread_num()
#endif


!$OMP DO SCHEDULE(STATIC,10) REDUCTION(+:force_long,v_sol_sol)     
    do i_part = 1,n_part  !n_mon_tot= brushes + droplet/melt
          i_dummy = ff_list(0,i_part)
          i_type = a_type(i_part)
# elif BIN_TYPE == 2  /* cell_list.f90 */
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(r_2_min, force_long, v_sol_sol, part_init_d, n_mon_tot, a_type, n_nei_cells, cell_of_part, lpart_in_cell, l_nei, l_eps, cell_neigh_ls, r0, r_cut2, inv_boundary, boundary, sigma_2_ss, e_shift, sig_long, mass, friction, v, r_cut_ss, inv_r_cut_ss, inv_range_2, r_cut_dpd_2)

#ifdef _OPENMP
    ith=omp_get_thread_num()
#endif

!$OMP DO SCHEDULE(STATIC,10) REDUCTION(+:force_long,v_sol_sol)

do i_part = part_init_d + 1, n_mon_tot ! loop over melt  
    i_cell = cell_of_part(i_part) ! get the cell of the particle i_part
    i_type = a_type(i_part)
# endif

          f_ipart(:) =0.
          q_part=i_part !dummy variable for OpenMP paralelization


# if BIN_TYPE == 0 || BIN_TYPE == 1
             do i_neigh = 1, i_dummy
                 j_part = ff_list(i_neigh,i_part)
                 j_type = a_type(j_part)



# elif BIN_TYPE == 2  /* cell_list.f90 */
           do j_dummy = 1, n_nei_cells ! loop neighbor cells
                j_cell = cell_neigh_ls(i_cell,j_dummy) ! j_cell = neighbour cell
        
                j_part = lpart_in_cell(j_cell)
                !DEBUG
                !print*,"i_cell,j_cell",i_cell,j_cell
                !/DEBUG
                do while(j_part .ne. 0) ! loop over neighbour particles in j_cell
    
                !DEBUG
                !print*,"i_part,j_part",i_part,j_part

                !DEBUG
                if(i_part .eq. j_part) print*,"ERROR: i_part,j_part",i_part,j_part
                j_type = a_type(j_part)
                     if( j_type .ne. 3 ) exit ! If particle not solvent, go to next cell. Cell linked list is done,
                    ! such that the last  particles in a cell are always solvent. If this particle is not solvent,
                   ! there is no more solvent in this cell. Go to next cell                                 
# endif
             if (j_part .le. part_init_d ) print*,"ERROR j_part in brush",j_part           
                  
              ! HPC

              delta_r(1) = r0(1,i_part) - r0(1,j_part)
              delta_r(2) = r0(2,i_part) - r0(2,j_part)
              delta_r(3) = r0(3,i_part) - r0(3,j_part)

              !----- PBC ----
              ! HPC

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              if(r_2 .le. 0.000001) then
                  r_2=0.01
                  print*,"Error, distance between particles is 0"
                  print*,"Particles",i_part,j_part
                  print*,r0(:,i_part)
                  print*,r0(:,j_part)
              end if
              !-----  check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then
#           ifdef FORCE_SWITCH_ON 
                  r_2 = max(r_2,r_2_min) 
#           endif                 

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                 ! IF SOL SOL INT ARE SOFT
#if SOL_SOL_INT == 2
                  r_dummy =  (inv_sqrt_r_2 - inv_r_cut_ss )
                v_sol_sol = v_sol_sol + 0.5 * l_eps * r_cut_ss * r_2 *r_dummy**2 
                r_dummy = - l_eps * r_dummy 
                ! ELSE IF SOL SOL INTERACTIONS ARE LJ: 
#elif SOL_SOL_INT == 1
                  r_61= sigma_2_ss*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  v_sol_sol = v_sol_sol + l_eps*pot_loc

                  r_dummy = l_eps*(-12*r_12+6*r_6)*inv_r_2
#endif /* SOL_SOL_INT*/
                  !! END IF

                  force_loc(1) = r_dummy*delta_r(1)
                  force_loc(2) = r_dummy*delta_r(2)
                  force_loc(3) = r_dummy*delta_r(3)
                  ! HPC

!DEBUG
if ( dot_product(force_loc,force_loc) .gt. 10.) then
    print*,"Problem in conservative force in calc_solv_solv"
    print*,"i_part,j_part",i_part,j_part
    print*,"r0(:,i_part), r0(:,j_part)",r0(:,i_part), r0(:,j_part)
    print*,"force",force_loc(:)
end if
!DEBUG

#       if BIN_TYPE == 0 ||  BIN_TYPE == 2/* binning.f90 */
                  f_ipart(1) = f_ipart(1) -  force_loc(1)
                  f_ipart(2) = f_ipart(2) -  force_loc(2)
                  f_ipart(3) = f_ipart(3) -  force_loc(3)
#       endif

                  force_long(1,j_part) = force_long(1,j_part) + force_loc(1)
                  force_long(2,j_part) = force_long(2,j_part) + force_loc(2)
                  force_long(3,j_part) = force_long(3,j_part) + force_loc(3)

#   if SYMMETRY == 1 /* bulk */                  
! Does not work with OMP
! Potential part of press tensor (see viscosity.f90) ~ virial contribution
         do i = 1,3
             do j = 1,3               
                 press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)              
             end do 
         end do 
!WRONG?  end if
    
#   endif

!
!  --- DPD calculation if DPD has not its own cut-off  
!

#if THERMOSTAT == 0
#     ifndef DPD_CUT_OFF
#           if SYMMETRY == 0
#                  if PINNED == 1
                  if(i_type.eq.3) cycle       !  Excluding  fixed particles
                  if(i_type.eq.4) cycle       ! 
                  if(j_type.eq.3) cycle       ! 
                  if(j_type.eq.4) cycle       ! 
#                  endif 
#         endif /* SYMMETRY */

                  call dpd_forces(inv_sqrt_r_2,force_long, sig_long)

#       endif /* not DPD_CUT_OFF */
#endif /*THERMOSTAT = 0 */


#       if SYSTEM == 2  || SYSTEM  == 3 /* charged systems */

! Coulomb potential and interaction in REAL space 

                  call ewald_real()

#       endif

              end if ! if the particle is inside the interaction sphere 

!
! DPD has its own cutoff radius 
!

#ifdef DPD_CUT_OFF
#       if THERMOSTAT == 0
#                      if SYMMETRY == 0
#                              if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3
              
#                                  if PINNED == 1

              if(i_type.eq.3) cycle       !  Excluding  fixed particles
              if(i_type.eq.4) cycle       ! 
              if(j_type.eq.3) cycle       ! 
              if(j_type.eq.4) cycle       ! 
#                                  endif 

#                              endif
#                       endif 

              if(r_2 < r_cut_dpd_2 ) then 
                  
!DEBUG
!print*,"Entro a dpd_forces"
                  call dpd_forces(inv_sqrt_r_2,force_long, sig_long)
!print*,"SAlgo dpd_forces"
              end if

#       endif /* THEMOSTAT = 0 */
#endif /* DPD_CUT_OFF */

#if BIN_TYPE == 0 || BIN_TYPE == 1
          end do ! loop over particles neighbors (j_part)
#elif BIN_TYPE == 2
                   
                    j_part = l_nei(j_part)

                end do ! loop over particles in j_cell
            end do ! loop over neighbor cells 

           j_part =  lpart_in_cell(i_cell)
            do while(j_part .ne. 0)  ! loop over particles in the same cell

             j_type = a_type(j_part)
             if( j_type .ne. 3) exit

            if (i_part .lt. j_part) then !Count interactions just once
                            ! HPC

              delta_r(1) = r0(1,i_part) - r0(1,j_part)
              delta_r(2) = r0(2,i_part) - r0(2,j_part)
              delta_r(3) = r0(3,i_part) - r0(3,j_part)

              !----- PBC ----
              ! HPC

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then
#           ifdef FORCE_SWITCH_ON 
                  r_2 = max(r_2,r_2_min) 
#           endif                 

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                 ! IF SOL SOL INT ARE SOFT
#if SOL_SOL_INT == 2
                  r_dummy =  (inv_sqrt_r_2 - inv_r_cut_ss )
                v_sol_sol = v_sol_sol + 0.5 * l_eps * r_cut_ss * r_2 *r_dummy**2 
                r_dummy = - l_eps * r_dummy 
                ! ELSE IF SOL SOL INTERACTIONS ARE LJ: 
#elif SOL_SOL_INT == 1
                  r_61= sigma_2_ss*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  v_sol_sol = v_sol_sol + l_eps*pot_loc

                  r_dummy = l_eps*(-12*r_12+6*r_6)*inv_r_2
#endif /* SOL_SOL_INT*/
                  !! END IF

                  force_loc(1) = r_dummy*delta_r(1)
                  force_loc(2) = r_dummy*delta_r(2)
                  force_loc(3) = r_dummy*delta_r(3)

                  !DEBUG
if ( dot_product(force_loc,force_loc) .gt. 10.) then
    print*,"Problem in conservative force in calc_solv_solv"
    print*,"i_part,j_part",i_part,j_part
    print*,"r0(:,i_part), r0(:,j_part)",r0(:,i_part), r0(:,j_part)
    print*,"force",force_loc(:)
end if
!DEBUG


                  f_ipart(1) = f_ipart(1) -  force_loc(1)
                  f_ipart(2) = f_ipart(2) -  force_loc(2)
                  f_ipart(3) = f_ipart(3) -  force_loc(3)

                  force_long(1,j_part) = force_long(1,j_part) + force_loc(1)
                  force_long(2,j_part) = force_long(2,j_part) + force_loc(2)
                  force_long(3,j_part) = force_long(3,j_part) + force_loc(3)

#   if SYMMETRY == 1 /* bulk */                  
! Does not work with OMP
! Potential part of press tensor (see viscosity.f90) ~ virial contribution
         do i = 1,3
             do j = 1,3               
                 press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)              
             end do 
         end do 
!WRONG?  end if
    
#   endif
 
!
!  --- DPD calculation if DPD has not its own cut-off  
!

#if THERMOSTAT == 0
#     ifndef DPD_CUT_OFF
#           if SYMMETRY == 0
#                 if PINNED == 1
                  if(i_type.eq.3) cycle       !  Excluding  fixed particles
                  if(i_type.eq.4) cycle       ! 
                  if(j_type.eq.3) cycle       ! 
                  if(j_type.eq.4) cycle       ! 
#                 endif 
#         endif /* SYMMETRY */

                  call dpd_forces(inv_sqrt_r_2,force_long, sig_long)

#       endif /* not DPD_CUT_OFF */
#endif /*THERMOSTAT = 0 */


#       if SYSTEM == 2  || SYSTEM  == 3 /* charged systems */

! Coulomb potential and interaction in REAL space 

                  call ewald_real()

#       endif

              end if ! if the particle is inside the interaction sphere 

!
! DPD has its own cutoff radius 
!

#ifdef DPD_CUT_OFF
#       if THERMOSTAT == 0
#                      if SYMMETRY == 0
#                              if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3
              !                               Exclude heads from DPD calculation 
#                                  if PINNED == 1

              if(i_type.eq.3) cycle       !  Excluding  fixed particles
              if(i_type.eq.4) cycle       ! 
              if(j_type.eq.3) cycle       ! 
              if(j_type.eq.4) cycle       ! 
#                                  endif 

#                              endif
#                       endif 

              if(r_2 < r_cut_dpd_2 ) then 
                  
!print*,"Entro a dpd_forces"
                  
                  call dpd_forces(inv_sqrt_r_2,force_long, sig_long)
!print*,"Salgo a dpd_forces"

              end if

#       endif /* THEMOSTAT = 0 */
#endif /* DPD_CUT_OFF */

           
                
                end if !  i_part < j_part  
            j_part = l_nei(j_part) !gets next particle in this cell     
            end do  ! loop over neighbor particles in i_cell, j_part
 
#endif


#       if BIN_TYPE == 0 || BIN_TYPE == 2/* binning.f90 */
 
          force_long(1,i_part) = force_long(1,i_part)  + f_ipart(1)
          force_long(2,i_part) = force_long(2,i_part)  + f_ipart(2)
          force_long(3,i_part) = force_long(3,i_part)  + f_ipart(3)

#       endif

#   if THERMOSTAT == 1 /* Langevin */
! COMMENT to check NVE ensemble          
 call lgv_forces(force_long,sig_long)
#   endif /* Langevin */

#if BIN_TYPE == 0 || BIN_TYPE == 1
end do ! loop over particles

!$OMP END DO 
!$OMP END PARALLEL 
!END PARALLEL ZONE

#elif BIN_TYPE == 2

      end do ! loop over all melt particles

!$OMP END DO 
!$OMP END PARALLEL 
!END PARALLEL ZONE


#endif

!DEBUG
!DEBUG

!do i_part = part_init_d + 1, n_mon_tot
!    print*,force_long(:,i_part)
!end do
        !print*,"v_sol_sol",v_sol_sol
!DEBUG
!DEBUG


#ifdef FLUKT 
    if((n_mon_d.or.n_chain_d.or.n_mon_e.or.n_chain_e).gt.2) then
    write(*,*) "the file fort.555 gonna be too big! undefine flag FLUKT or decrease number of particles4 or -3"
    stop
    endif
    if (mod(i_time, 100).eq.0) then
    write(555,'(1i,18f15.4)') i_time,  r0(part_init_d+1:n_mon_tot, :), & 
                              v(:,part_init_d+1:n_mon_tot), force_long(part_init_d+1:n_mon_tot, :)
    end if
#endif
#                 if BIN_TYPE == 1                  
                  v_sol_sol = 0.5*v_sol_sol
#                 endif 
#if SYSTEM == 2 || SYSTEM == 3
        v_coul = v_coul + v_coul_self
#endif

! Debugging
if(debug) then
!print '(3f15.5)', force(:,1:n_chain*n_mon:n_mon)
          print *,"v_sol_sol=",i_time,v_sol_sol/dble(n_mon_tot)
r_dummy =     sum( sum(force(:,:)**2,dim=1) ,dim =1 )
print *,"[fluid_fluid]quad_m_force", i_time,sqrt(r_dummy)/dble(n_mon_tot)
print *, "[fluid_fluid]V=",v_sol_sol
!do i_part=1,n_mon_tot
!write(88,'(3i,3f15.4)') i_time, i_part,a_type(i_part),force(:,i_part)
!end do
end if

!DEBUG
!print*,"Salgo a calc_solv_solv_force"
!print*,""
!!/DEBUG         !

end subroutine calc_solv_solv_force
