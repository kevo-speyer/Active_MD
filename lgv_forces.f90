
!! lgv_forces !! 
! * Computes random and friction forces for langevin thermostat 
! * It is assumed that the routine is called inside the fluid_fluid routine

subroutine lgv_forces(force_tmp)
    use commons
#ifdef _OPENMP
      use omp_lib
      use Par_Zig_mod, only: par_rnor,par_uni
#else
      use ziggurat, only: rnor,uni
#endif
#include "control_simulation.h" 
    implicit none
!    real(kind=8) , intent(in) :: inv_sqrt_r_2
    integer :: i
!    integer, intent(in) :: q_part
    real(kind=8), intent(inout) :: force_tmp(3,n_part)
    real(kind=8) :: g_rand,r_versor(3),delta_v(3),rrc
!    real(kind=8),intent(out) :: f_lang(3) 

!    *** Uniform distribution, claimed to be as good as the gaussian
! (fullfills 1st and second momentum properties of the random variable)
! See: Paul & Duenweg 

!       g_rand = sqrt(3.)*(2.*uni() - 1.) 


! ----  Random force computation:      

!print *, sig ; stop
!    ith=omp_get_thread_num()

    do i=1,3
#ifdef _OPENMP
       g_rand=par_rnor(ith)
!print*, jth,g_rand
#else
       g_rand=rnor()
#endif
       vec_dummy(i) = sig*g_rand*sqrt(mass(q_part))    
    end do
!print*,"q_part lgv",q_part
!original:      vec_dummy(1) = sig*rnor()*sqrt(mass(i_part))
!original:      vec_dummy(2) = sig*rnor()*sqrt(mass(i_part))
!original:      vec_dummy(3) = sig*rnor()*sqrt(mass(i_part))
   
!      force(:,i_part) =  force(:,i_part) + vec_dummy(:)*sqrt(mass(i_part))   

     
! ---- Dissipative force computation       
      vec_dummy(:) = vec_dummy(:) -friction(1)*mass(q_part)*v(:,q_part)
!      !ORIGINAL
#   ifdef LGVX_0
            vec_dummy(1) = 0. ! No thermostat in X
#   endif
     force_tmp(:,q_part) =  force_tmp(:,q_part) + vec_dummy(:)   
end subroutine lgv_forces 
