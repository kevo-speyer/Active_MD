subroutine new_dpd_fd()
! 
! Recalculation of Fd  after update of velocities. Only Fd from DPD thermostat is recalculated
!  * Also  
!
#include 'control_simulation.h'
use commons
!use ziggurat, only: rnor,uni
!use util ! debug
  implicit none
  real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc ! w_d and w_r now in mfa_commons.f90
  real (kind=8) :: T_inst
  integer (kind=8) :: n_trial
  logical, parameter :: debug=.false.
  
! ----- Substract first  the  current Fd to the total force 

        force(:,:)  = force(:,:) - force_d(:,:)

! ----   Init. Fd.

        force_d(:,:) = 0.

!
! New dissipative force calculation: 
! (using linked list calculated in binning.f90)

        do i_part = 1,n_mon_tot  !n_mon_tot= brushes + droplet/melt
            i_type =a_type(i_part)
            if(i_type.eq.1) cycle          ! Heads are excluded from this pair calculation !!!

            i_dummy = ff_list(0,i_part)
            do i_neigh = 1,i_dummy
                j_part = ff_list(i_neigh,i_part)
                j_type =a_type(j_part)
                if(j_type.eq.1) cycle       ! Excluding heads !!!

                ! -- Boundary conditions in x and y 

                delta_r(:) = r0(:,i_part) - r0(:,j_part)

                ! In the old good way ...

                do i_dim = 1,n_dim-1
                    delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)*inv_boundary(i_dim))
                end do
                r_2 = 0.
                do i_dim = 1,n_dim
                    r_2 = r_2 + delta_r(i_dim)**2
                end do


                ! -- If interaction takes place 

                if(r_2.lt.range_2(i_type,j_type)) then

                    !    *** Gaussian random distribution. sig^2 = 1 , <x> = 0

                    ! Here the complete story of different weights is reproduced 

                    ! --- Weight functions 

#if DPD_WEIGHT == 0
                    w_d = (1. - sqrt(r_2/range_2(i_type,j_type)))**2  
#endif      

!!note:  #if DPD_WEIGHT == 1 as the functions are constants, they are defined in init_params.f90
!!               and w_d=w_r for ever in the rest of the program
!!note:        w_d = 1.
!!note:  #endif      

#if DPD_WEIGHT == 2
                    w_d = sqrt( 1. - sqrt(r_2/range_2(i_type,j_type)) )  
#endif      
#if DPD_WEIGHT == 3
                    rrc= sqrt(r_2/range_2(i_type,j_type))
                    w_d = ( 1. - rrc )**0.25  
#endif      
#if DPD_WEIGHT == 4
                    rrc= sqrt(r_2/range_2(i_type,j_type))
                    w_d = ( ( 1. - rrc )**0.25 ) *100.* exp (-(1-rrc)**2)  
#endif      


                    r_versor(:) = delta_r(:)/sqrt(r_2)

                    !not in dpd_fd      w_r = sqrt(w_d)

                    delta_v(:) = (v(:,i_part) - v(:,j_part) )
                    r_dummy = delta_v(1)*r_versor(1) + delta_v(2)*r_versor(2) + delta_v(3)*r_versor(3)
                    vec_dummy(:) = -1.*friction(1)*w_d*r_dummy*r_versor(:)

                    force_d(:,i_part) =  force_d(:,i_part) + vec_dummy(:)   
                    force_d(:,j_part) =  force_d(:,j_part) - vec_dummy(:)   

                end if

            end do ! particles in the cell
        end do    ! particles

        !  --- Add the new Fd to the force 

        force(:,:) = force(:,:) + force_d(:,:) 

! Note:  For dpd thing we calculate the temperature only here (and not in verlet_velocities_)

!
! What follows has been copy/pasted from verlet_velocities: Es,T and the new accellerations are calculated only
! after the correction of the of the Dissipative force.  
!

do i_part = 1 , n_mon_tot
!    inv_mass = 1./mass(i_part)
 ! here the  velocity is updated until half the interval   
      v(:,i_part) = v(:,i_part) +  0.5*dt*force(:,i_part)*inv_mass(i_part) 
end do

! Update accelerations with the new force values 

    do i_part = 1,n_mon_tot
        a(1,i_part) = force(1,i_part)*inv_mass(i_part)
        a(2,i_part) = force(2,i_part)*inv_mass(i_part)
        a(3,i_part) = force(3,i_part)*inv_mass(i_part)
    end do

!
! rx to 'wall velocity' for brush heads, if there are brushes 
!
#    if WALL == 1 /* explicit wall */
    if(n_chain>0) then
        do i_chain = 1 ,n_chain
            if(i_chain<=n_chain/2) then
                v(1:2,1+n_mon*(i_chain-1)) = va_spring_twall(1:2)
            else  
                v(1:2,1+n_mon*(i_chain-1)) = -va_spring_twall(1:2)
            end if
            v(3,1+n_mon*(i_chain-1)) = 0.
        end do
    end if
#   endif



!(from corrector.f90)
!---  calculate kinetic energies
!
      t_fluid = 0.
      t_wall = 0.
      T_inst = 0.
!
    do i_dim= 1,n_dim
        do i_part = 1,n_mon_tot
            t_fluid = t_fluid + mass(i_part)*v(i_dim,i_part)**2

            ! Instantaneous temperature

            T_inst = T_inst + mass(i_part)*v(i_dim,i_part)**2
        end do

        !if(f_explicit_wall) then       
#   if WALL == 1 /* explicit wall */           
        do i_part = n_mon_tot+1,n_part
            t_wall = t_wall + mass(i_part)*v(i_dim,i_part)**2
        end do
#   endif
        !end if
    end do


! *** Write out Instantanous T and momentum  ! then directly yo a file WARN
!       if(mod(i_time,500).eq.0) then
!            T_inst_mean  = T_inst_mean +  T_inst / ((3.*real(n_mon_tot-n_chain) - 2.)*dt_2) 
!            print *,"T_inst=",i_time,T_inst
!            vec_dummy = 0.
!            do i_part = 1,n_part
!            vec_dummy(:) = vec_dummy(:) + mass(i_part)*dt*v(:,i_part)
!            end do
!            print '(a20,3f15.5)',"Total momentum P=", vec_dummy(:)/dt 
!       end if
!
      t_wall = t_wall/2.
      t_fluid = t_fluid/2.
!
      v_total = v_wall_wall + v_fluid_wall + v_fluid_fluid + v_intra_molec
!
      if(f_twall(n_dim).eq.2) then
       v_total=v_total-r0_twall(n_dim)*va_spring_twall(n_dim)/2
      end if
!
      t_total = t_wall + t_fluid
      e_total = v_total + t_total
!      
!  Temperature and energy calculation and writing out 
!

   time_ave_count_2 = time_ave_count_2 + 1
   if(time_ave_count_2.eq.10) then
       time_ave_count_2 = 0
       do i_dim = 1,n_dim
        r0_twall_1(i_dim) = r0_twall_1(i_dim) / n_time_ave
        f0_spring_1(i_dim) = f0_spring_1(i_dim)/n_time_ave
       end do
    if(mod(i_time,1).eq.0) then ! How frequently write out 
        write(60,'(i11,3f17.5)') i_time,v_total/real(n_part-n_chain),            &
        e_total/dble(n_part-n_chain),t_total/dble(n_part-n_chain)
       end if
       do i_dim = 1,n_dim
        r0_twall_1(i_dim) = 0.
        f0_spring_1(i_dim) = 0.
       end do
    end if


end subroutine new_dpd_fd

