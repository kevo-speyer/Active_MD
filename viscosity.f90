subroutine viscosity()
! Pressure tensor and all physical magnitudes involved in viscosity calculations
! are calculated here.
!   
use commons
implicit none
logical, parameter :: debug=.false.
integer :: i,j
real(kind=8) :: delta_f(3), visc_t = 0., norm_press_t_corr(3,3) = 0.
real(kind=8),save :: press_tensor_0(3,3),press_ten_f(3,3)
 
#if SYMMETRY == 1 /* Bulk simulation */

! Add Kinetic part of pressure tensor

press_tensor_0(:,:) = 0.


press_ten_f(:,:) =  press_tensor(:,:)

    do i_part= 1,n_mon_tot
        do i = 1,3
            do j = 1,3
                press_tensor(i,j) =  press_tensor(i,j) + mass(i_part)*v(i,i_part)*v(j,i_part )
                press_tensor_0(i,j) =  press_tensor_0(i,j) + mass(i_part)*v(i,i_part)*v(j,i_part ) !deb
            end do 
        end do 
    end do 

! In the current definition, normalize also with the volume of the system


      press_tensor(:,:) = press_tensor(:,:) * inv_V
      ! press_tensor_0 :: Kinetic part of compress. factor.
      press_tensor_0(:,:) = press_tensor_0(:,:) !  /(temp*dble(n_part)) !deb

      
    
! [NOW in each routine which compute forces !] Potential part of the pressure tensor     

!!! ************************* LOOP OVER PAIRS WITH linked cell list calculated in binning.f90    
!! do i_part = 1,n_mon_tot  !n_mon_tot= brushes + droplet/melt
!!       i_type =a_type(i_part)
!!
!!       i_dummy = ff_list(i_part,0)
!!       do i_neigh = 1,i_dummy
!!        j_part = ff_list(i_part,i_neigh)
!!       j_type =a_type(j_part)
!!
!!! -- Boundary conditions in x and y  AND Z
!!
!!     delta_r(:) = r0(i_part,:) - r0(j_part,:)
!!
!!        do i_dim = 1,n_dim
!!         delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)/boundary(i_dim))
!!        end do
!!     
!!        delta_f(:) = force(i_part,:) - force(j_part,:)
!!
!!        do i = 1,3
!!        do j = i+1,3
!!        
!!        press_tensor(i,j) =  press_tensor(i,j) + delta_f(i)*delta_r(j)
!!        press_tensor(j,i) =  press_tensor(j,i) + delta_f(j)*delta_r(i)
!!
!!        end do 
!!        end do 
!!
!!     end do ! particles in the cell
!!  end do    ! particles
!!  
!!  if (i_time==1) then
!!
!!      press_tensor_0(:,:) =  press_tensor(:,:)
!!
!!        do i = 1,3
!!        do j = i+1,3
!!        
!!        norm_press_t_corr(i,j) =  press_tensor_0(i,j)*press_tensor_0(i,j)
!!        norm_press_t_corr(j,i) =  press_tensor_0(j,i)*press_tensor_0(j,i)
!!
!!        end do 
!!        end do 
!!
!!      
!!  end if

! Write to a file each component of pressure tensor vs. T and the <>: sum/6  
! (from the integral of the correlation function of the press tensor, shear viscosity can be calculated)  

 if (mod(i_time,10) == 0 ) then    
  write(unit=76,fmt='(f18.4,7e18.5)') r_time,                                                 &  
     press_tensor(1,2),  press_tensor(1,3),   press_tensor(2,3),   press_tensor(2,1),   & !components
       press_tensor(3,1),   press_tensor(3,2) 

  write(unit=77,fmt='(f18.4,3g18.5)') r_time,                                                 &  
     press_tensor(1,1),  press_tensor(2,2),   press_tensor(3,3)

  write(unit=78,fmt='(f18.4,6g18.5)') r_time,                                                 &  
     press_tensor_0(1,1),  press_tensor_0(2,2),   press_tensor_0(3,3), &
     press_ten_f(1,1),  press_ten_f(2,2), press_ten_f(3,3)
   end if
       


! Compute mean value of pressure tensor

      press_tensor_mean(:,:) =  press_tensor_mean(:,:) +  press_tensor(:,:) 

       
#endif /*Bulk: SYMMETRY == 1 */
! Note: the pressure tensor is calculated WITHOUT the factor 1/V, because is then cancelled out
! in the correlation function for viscosity calculation

! The viscosity is calculated as (V/kT) * integral < press_tensor(0)*press_tensor(t) > dt 
! for each non-linear element of pressure tensor. 
!
! WARN: correlation function of stress tensor is NEEDED to do that.

! visc_t =  visc_t + ( press_tensor_0(1,2)*press_tensor(1,2) /norm_press_t_corr(1,2)+   &  
!                      press_tensor_0(1,3)*press_tensor(1,3) /norm_press_t_corr(1,3)+   &
!                      press_tensor_0(2,3)*press_tensor(2,3) /norm_press_t_corr(2,3)+   &
!                      press_tensor_0(2,1)*press_tensor(2,1) /norm_press_t_corr(2,1)+   &  
!                      press_tensor_0(3,1)*press_tensor(3,1) /norm_press_t_corr(3,1)+   &
!                      press_tensor_0(3,2)*press_tensor(3,2) /norm_press_t_corr(3,2)    &
!                      ) /6.*dt/temp
! check
!print*, 'visc', r_time , visc_t 
                      
end subroutine viscosity 

