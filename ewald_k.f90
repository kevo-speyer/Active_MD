    subroutine ewald_k(mode)
! This routine computes, ewald sum implementation of the Coulomb interactions
! in the real space

#include 'control_simulation.h'
   use commons
   implicit none 
   real(kind=8) :: u_cou_k,rk_2
   real(kind=8), allocatable, save :: k_vecs_2(:),inv_k_vecs_2(:),cphase(:,:)
   real(kind=8),save :: delta_k(3),fac_alpha,inv_alpha_c,tot_q,k_c
   real(kind=8), allocatable, save :: re_rho_k(:),im_rho_k(:),force_c_k(:,:),k_vecs(:,:)
   real(kind=8) :: f_k(3),lphase
   integer ,intent(in) :: mode
   integer :: nkx,nky,nkz,k_index
   integer, save :: n_kmax_d,Nk,n_kmax_z,n_kmax_y,n_kmax_x
!   real(kind=8), save :: u_cou_k_head
   logical , parameter :: debug=.false.

#if SYSTEM == 2  || SYSTEM == 3
select case (mode)

case(0)  ! Init k variables 

print '(/a/) ' , "  *  Initializing variables for Ewald in k-space " 

! NOTE: We calculate all the K's from: 0,kz_max ; -kx_max,kx_max, -ky_max,k_ymax,
! Excluding k = (0,0,0)  and other points 
! The plane z=0 has fewer points than all other planes in the positive kz part
! of K-space
!  [total # of ks] =       [  kz=0 plane ]    +  [ kz >0 planes ]      
! [ kz >0 planes ] : #k = n_kmax*(2*n_kmax+1)^2
! 
![ kz=0 plane ]  : More complicated. The lower triangle must be excluded and
! also : half of the points lying in the rect y = -x . In the current
! implemmentation, the points with positive ky are omitted. 

#   if SYMMETRY == 1 /* if bulk geometry */
        delta_k(:)  = 2.*pi*inv_boundary(:)
        n_kmax_z = n_kmax
        z_slab = 0. 
#   elif SYMMETRY == 0  /* if we have a channel or a slab */
        z_slab = 2. 

        delta_k(1:2)  = 2.*pi*inv_boundary(1:2)
        delta_k(3) = 2.*pi/( boundary(3) + z_slab*boundary(3) ) 
#   endif
        k_c = (delta_k(1)+delta_k(2))* dble(n_kmax)/2.

        n_kmax_x = int( k_c/delta_k(1) )
        n_kmax_y = int( k_c/delta_k(2) )
        n_kmax_z = int( k_c/delta_k(3) )

!        n_kmax_z = n_kmax*int(z_slab + 1.) ! n_kmax*z_slab

        print '(/a,f8.3)' ,     "  *  K_max = " ,k_c
        print '(a,3i4)' ,       "  *  (nx,ny,nz)_max = " ,n_kmax_x,n_kmax_y,n_kmax_z
        print '(a,f8.3,a)' ,    "  *  Taking a slab separation = " , z_slab , " Lz"
        print '(a,i6)' ,        "  *  Points in Z = " ,n_kmax_z
        print '(a,3f12.5)' ,    "  *  delta_k=  " , delta_k(:)
        print '(a/)' ,          "  *  No PBC in Z for Ewald in k-space" 

        k_index = 0
        do nkz = 0,n_kmax_z*3
            do nky = -n_kmax_y,n_kmax_y
                do nkx = -n_kmax_x,n_kmax_x
                    rk_2 = (dble(nkx)*delta_k(1) )**2 + (dble(nky)*delta_k(2) )**2 + (dble(nkz)*delta_k(3) )**2
                    if (rk_2 >  k_c**2) cycle ! spherical symmetry
                    if (nkz == 0 .and. nky == 0 .and. nkx == 0 )  cycle  ! exclude k=0
                    if (nkz == 0 .and. nkx+nky < 0 ) cycle   ! exclude lower triangular part from the points in the plane kz =0
                    if (nkz == 0 .and. -nkx==nky .and. nky > 0   ) cycle   ! exclude half of the diagonal from the points  
                    k_index = k_index + 1
                    !                print*,"k=",k_index,nkx,nky,nkz
                    print '(a,i6,3f12.5)',"k=",k_index,dble(nkx)*delta_k(1),dble(nky)*delta_k(2),dble(nkz)*delta_k(3)
                end do
            end do
        end do

        Nk = k_index

    print '(/a,i6/) ' , "  *  K-vectors to be used for Ewald = " ,Nk

        allocate ( re_rho_k(Nk) ,            &
                   im_rho_k(Nk ) ,           &
                   force_c_k(3,n_mon_tot)  , & 
                    k_vecs(Nk,3) )
 
        allocate(  k_vecs_2(Nk),inv_k_vecs_2(Nk),cphase(n_mon_tot,Nk) )

        fac_alpha = 1./(4.*alpha_c)
 

        inv_V=inv_boundary(1)*inv_boundary(2)/( boundary(3)*(1.+z_slab) )

        inv_alpha_c = 1./alpha_c

        tot_q = sum(q(:),dim=1)

        print '(/a,f12.3/)', "  *  Total charge of the system = ", tot_q

!
! ----- Calculation of  k_vecs to be used for the rest of the run
!

        k_index = 0
        do nkz = 0,n_kmax_z
            do nky = -n_kmax_y,n_kmax_y
                do nkx = -n_kmax_x,n_kmax_x

                    rk_2 = (dble(nkx)*delta_k(1) )**2 + (dble(nky)*delta_k(2) )**2 + (dble(nkz)*delta_k(3) )**2

                    if (rk_2 > k_c**2) cycle ! spherical symmetry
                    if (nkz == 0 .and. nky == 0 .and. nkx == 0 )  cycle  ! k=0
                    if (nkz == 0 .and. nkx+nky < 0 ) cycle   ! exclude lower triangular part from the points in the plane kz =0
                    if (nkz == 0 .and. -nkx==nky .and. nky > 0   ) cycle   ! exclude half of the diagonal from the points  

                    k_index = k_index + 1
                    k_vecs(k_index,:) =  delta_k(:)*dble((/ nkx,nky,nkz /))                
                    k_vecs_2(k_index) = dot_product(k_vecs(k_index,:),k_vecs(k_index,:)) 
                    inv_k_vecs_2(k_index) = 1./k_vecs_2(k_index)


                end do
            end do
        end do



case(1)  ! K-space calculation

! ------ Rho-k: Calculation of Re_rho_k and Im_rho_k for positive k's

         u_cou_k = 0.
         re_rho_k(:) = 0.
         im_rho_k(:) = 0.

         do i_part = 1,n_mon_tot
             do k_index = 1 , Nk

                 lphase = dot_product( k_vecs(k_index,:),r0(:,i_part) )
                 cphase(i_part,k_index)= lphase

                 re_rho_k(k_index) = re_rho_k(k_index) + q(i_part)*cos(lphase)
                 im_rho_k(k_index) = im_rho_k(k_index) + q(i_part)*sin(lphase)

             end do  
         end do


! Ucoul_k calculation

         do k_index = 1 , Nk
!                k_vecs_2 = dot_product (k_vecs(k_index,:),k_vecs(k_index,:) )
!                inv_k_vecs_2 = 1./k_vecs_2

                u_cou_k = u_cou_k + inv_k_vecs_2(k_index)*exp(-k_vecs_2(k_index)*fac_alpha)*   &
                (re_rho_k(k_index)**2 + im_rho_k(k_index)**2 )
        end do

                u_cou_k = u_cou_k*4.*pi*inv_V

! Note:  (re_rho_k*re_rho_k+im_rho_k*im_rho_k) == |rho(k)|^2
!
! Coulomb force calculation in k-space
!
!    force_c_k(:,:) = 0.


    do i_part = 1, n_mon_tot
        f_k(:) = 0. ! Accummulates the interaction for a given i_part
        do k_index = 1,Nk 

!                k_vecs_2 = dot_product (k_vecs(k_index,:),k_vecs(k_index,:) )
!                inv_k_vecs_2 = 1./k_vecs_2

!NOTE:  (re_rho_k*re_rho_k+im_rho_k*im_rho_k) == |rho(k)|^2

!                cphase= dot_product( k_vecs(k_index,:),r0(:,i_part) )

! Contribution of  k_s

            f_k(:) = f_k(:) + k_vecs(k_index,:)*inv_k_vecs_2(k_index)*exp(-k_vecs_2(k_index)*fac_alpha)*  &
                (  cos(cphase(i_part,k_index))*im_rho_k(k_index) - sin(cphase(i_part,k_index))*re_rho_k(k_index) )

        end do

! NOTE: As fi(-k) = fi(k). fi is symmteric in k, we add a factor x2 to account
! for the other half of the K-space considered. This factor is already added ommitting the original
! factor 1/2.
! Therefore, the calculation is actually done for +kz and the factor x2 accounts
! for all the points with -kz
!       This is also true for the plane z=0 where the calculation is done for
!       upper triangle of the plane.

        force_c_k(:,i_part) = -8.*pi*q(i_part)*f_k(:)*inv_V 
        force(:,i_part)  =  force(:,i_part) + force_c_k(:,i_part)  

    end do ! i_part


! Add to the total v_coul potential energy           

                v_coul = v_coul + u_cou_k

end select

        if (debug) then 
            print *,"[ewald_k]",i_time,u_cou_k
            print *,"[ewald_k]force",sum(force_c_k(:,:)**2*inv_N,dim=1)
        end if
#endif 
    end subroutine ewald_k
