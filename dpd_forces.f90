subroutine dpd_forces(inv_sqrt_r_2,force_tmp)
    use commons
#ifdef _OPENMP
      use omp_lib
      use Par_Zig_mod, only: par_rnor,par_uni
#else
use ziggurat, only: rnor,uni    
#endif
    
#include "control_simulation.h" 

    implicit none
    real(kind=8) , intent(in) :: inv_sqrt_r_2
    real(kind=8) :: g_rand,r_versor(3),delta_v(3),rrc
    real(kind=8) , intent(inout) :: force_tmp(3,n_part)

    integer :: i,j 
!    integer, intent(in) :: q_part
!    integer , intent(in)      :: i_type,j_type
!
!   -------- DPD forces calculation : Version inside the force caculation  
!
! WARN: check-out if works with different masses as is. 

!    *** Gaussian random distribution. sig^2 = 1 , <x> = 0

!      g_rand(:) = (/ rnor(),rnor(),rnor() /)    

#ifdef _OPENMP
          g_rand = par_rnor(ith)
#else
          g_rand =  rnor()      
#endif
         
!    *** Uniform distribution, claimed to be as good as the gaussian
! (fullfills 1st and second momentum properties of the random variable)

!       g_rand = sqrt(3.)*(2.*uni() - 1.) 

!       print *,"random",g_rand(:)
!      if(debug)      print*,"random",g_rand 

! --- Weight functions 

#if DPD_WEIGHT == 0
#       ifndef DPD_CUT_OFF
             w_d = (1. - sqrt(r_2*inv_range_2(i_type,j_type)))  
#       else      
             w_d = (1. - sqrt(r_2/r_cut_dpd_2) )
#       endif      
          w_r = w_d
          w_d = w_d*w_d
#endif      

!!note:  #if DPD_WEIGHT == 1 as the functions are constants, they are defined in init_params.f90
!!               and w_d=w_r for ever in the rest of the program
!!note:        w_d = 1.
!!note:  #endif      

#if DPD_WEIGHT == 2
      w_d = sqrt( 1. - sqrt(r_2/range_2(i_type,j_type)) )  
      w_r = sqrt(w_d)
#endif      
#if DPD_WEIGHT == 3
      rrc= sqrt(r_2/range_2(i_type,j_type))
      w_d = ( 1. - rrc )**0.25  
      w_r = sqrt(w_d)
#endif      
#if DPD_WEIGHT == 4
      rrc= sqrt(r_2/range_2(i_type,j_type))
      w_d = ( ( 1. - rrc )**0.25 ) *100.* exp (-(1-rrc)**2)  
      w_r = sqrt(w_d)
#endif      

! ----  Random force computation:      

          r_versor(:) = delta_r(:)*inv_sqrt_r_2
         vec_dummy(:) = sig*w_r*g_rand*r_versor(:) 
!        print*,"q,j:",q_part,j_part,sig

!print*,"r_v",r_versor
!print*,"v_d",vec_dummy

!#       if BIN_TYPE == 1             
!             vec_dummy(:) = 0.5*vec_dummy(:)
!#       endif


#   if BIN_TYPE == 0 
      force_tmp(:,q_part) =  force_tmp(:,q_part) + vec_dummy(:)   
      force_tmp(:,j_part) =  force_tmp(:,j_part) - vec_dummy(:)   
#   elif BIN_TYPE == 1
      if (q_part>j_part) then 
      force_tmp(:,q_part) =  force_tmp(:,q_part) + vec_dummy(:)   
      force_tmp(:,j_part) =  force_tmp(:,j_part) - vec_dummy(:)   
      endif
#   endif
     
! Dissipative force computation      

      delta_v(:) = v(:,q_part) - v(:,j_part)

      r_dummy = delta_v(1)*r_versor(1) + delta_v(2)*r_versor(2) + delta_v(3)*r_versor(3)

      vec_dummy(:) = -1.*friction(1)*w_d*r_dummy*r_versor(:)
        
!#       if BIN_TYPE == 1             
!             vec_dummy(:) = 0.5*vec_dummy(:)
!#       endif


      force_tmp(:,q_part) =  force_tmp(:,q_part) + vec_dummy(:)   
#     if BIN_TYPE == 0      
      force_tmp(:,j_part) =  force_tmp(:,j_part) - vec_dummy(:)   
#     endif      


#if SYMMETRY == 1
! potential part of press tensor (see viscosity.f90) ~ virial contribution
        do i = 1,3
            do j = 1,3
                
                press_tensor(i,j) =  press_tensor(i,j) + vec_dummy(i)*delta_r(j)
             
            end do 
        end do 
#endif

end subroutine dpd_forces
