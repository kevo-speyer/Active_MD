! Routine that introduces an external driving periodical force on the grafted chains
subroutine metronome(mode_metro)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_metro
real(kind=8) :: r_last(3)=0., r_last_mod=0., r_2rel(3)=0., f_ext=0., k_old=0., hard=0. 
real(kind=8) :: aux1(3)=0., aux2(3)=0., aux3(3)=0., d1(3)=0., d2(3)=0., cosine=0., sine=0., S=0., invd1d2=0., kappa=0.
!real(kind=8), ALLOCATABLE :: curv(:)
! , beta1=80., beta2=88.
real(kind=8), SAVE :: cos_lim1, cos_lim2, sin_lim
integer :: l,i 
integer, ALLOCATABLE, SAVE :: kick(:)
real(kind=8), ALLOCATABLE, SAVE :: curv(:)
integer, SAVE :: count_local 


#ifdef ACTIVE_BRUSH
select case (mode_metro)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with active polymer"
print*,""
ALLOCATE( cos_mem(n_chain) )
ALLOCATE( kick(n_chain) )
ALLOCATE( curv(n_chain) )
cos_mem = 0.
kick = 0
curv = 0.
count_local = 1
cos_th=0.
cos_lim1 = cos(beta1*pi/180.)
cos_lim2 = cos(beta2*pi/180.)
print*, cos_lim2, cos_lim1
sin_lim = sin(beta1*pi/180.)

if (beta1.gt.beta2) then
    print *, " **  WARNING: beta1 is greater than beta2, this will render active mode useless"
    print *, "(beta1=",beta1,") > (beta2=",beta2,")"
    print *, " ** Exiting now..."
    STOP
endif

!!! Case 1 & 2 were commented for easier debugging

!!!case(1)  ! Apply forces where angle condition requires
!!!! We'll compare 'cos(theta)' rather than 'theta'
!!!! It's biyective in [0,pi]->[-1,1], and cheaper to calclulate
!!!!cos_lim1 = abs( cos(pi/2.- beta*pi/180.) ) 
!!!f_ext = 10000. 
!!!
!!!
!!!do l = 1, n_chain
!!!    !r_last(:) = r0(:, l*n_mon)-r0(:, 1+(l-1)*n_mon)
!!!    !r_last_mod = sqrt(dot_product(r_last,r_last))
!!!    !cos_th = r_last(1)/r_last_mod        
!!!    r_2rel = r0(:, 3+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)         !REDO!!!
!!!    cos_th = r_2rel(1)/SQRT(DOT_PRODUCT(r_2rel,r_2rel))
!!!    if( (cos_mem(l).le.cos_lim1 ).and.(cos_th.ge.cos_lim1) ) then 
!!!        kick(l) = kick(l) + 15
!!!        if(l.le.n_chain/2) then  ! top wall 
!!!            force(1, 2+(l-1)*n_mon) = force(1, 2+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
!!!            force(1, 3+(l-1)*n_mon) = force(1, 3+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
!!!            force(3, 2+(l-1)*n_mon) = force(3, 2+(l-1)*n_mon) + f_ext*(cos_lim1)!*sign(1.,cos_th)
!!!            force(3, 3+(l-1)*n_mon) = force(3, 3+(l-1)*n_mon) + f_ext*(cos_lim1)!*sign(1.,cos_th)
!!!        else if (l.gt.n_chain/2) then !bottom wall
!!!            force(1, 2+(l-1)*n_mon) = force(1, 2+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
!!!            force(1, 3+(l-1)*n_mon) = force(1, 3+(l-1)*n_mon) + f_ext*(sin_lim)!*sign(1.,cos_th)
!!!            force(3, 2+(l-1)*n_mon) = force(3, 2+(l-1)*n_mon) - f_ext*(cos_lim1)!*sign(1.,cos_th)
!!!            force(3, 3+(l-1)*n_mon) = force(3, 3+(l-1)*n_mon) - f_ext*(cos_lim1)!*sign(1.,cos_th)
!!!        end if
!!!
!!!    end if
!!!
!!!    cos_mem(l) = cos_th  ! Store new memory for next step
!!!
!!!End do ! end do l


!!!case(2)
!!!! cos_lim2 < cos_lim1 for this to work
!!!
!!!!print *, cos_lim1, cos_lim2
!!!
!!!do l = 1, n_chain
!!!
!!!    r_2rel = r0(:, 3+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)  
!!!    cos_th = r_2rel(1)/SQRT(DOT_PRODUCT(r_2rel,r_2rel))
!!!    k_old = k0(1+(l-1)*n_mon)  ! storing old k
!!!
!!!    if( cos_th.gt.cos_lim1 ) then  ! the "more horizontal" (1) section
!!!        k0(1+(l-1)*n_mon) = 5*k_bend
!!!    else if( cos_th.lt.cos_lim2 ) then ! the "more vertical" (2) section 
!!!        k0(1+(l-1)*n_mon) = k_bend
!!!    else ! the middle section
!!!        if( cos_mem(l).le.cos_lim2 ) then ! if it came from 2
!!!            k0(1+(l-1)*n_mon) = k_bend  ! soft k
!!!        else if( cos_mem(l).gt.cos_lim1 ) then ! if it came from 1
!!!            k0(1+(l-1)*n_mon) = 5*k_bend   ! hard k
!!!        else ! if it came from the middle
!!!            k0(1+(l-1)*n_mon) = k_old   ! keep old k
!!!        end if
!!!    end if
!!!
!!!    cos_mem(l) = cos_th  !saving old position
!!!end do

case(3)
! cos_lim2 < cos_lim1 for this to work
! Curvature check, testing here, looking for a more permanent location...
do l = 1, n_chain
    kappa = 0.
    do i = 1, n_mon-2

        aux1 = r0(:,i+(l-1)*n_mon) !u
        aux2 = r0(:,i+1+(l-1)*n_mon) !v
        aux3 = r0(:,i+2+(l-1)*n_mon) !w
        d1 = aux2-aux1 !v-u
        d2 = aux3-aux2 !w-v
        invd1d2 = 1./(NORM2(d1)*NORM2(d2))
        cosine = DOT_PRODUCT(d1,d2)*invd1d2
        S = NORM2(aux1-aux3)
        sine = NORM2(cross(d1,d2))*invd1d2 !SQRT(1 - cosine**2)
        !print*, S, sine
        !if( sine<0.001) then
        !    curv(l) = curv(l) + 1000.
        !end if
        kappa = kappa + .5*S/sine
        !curv(l) = curv(l) + .5*S/sine
        !curv(l) = curv(l) + 2.*sine/S !+.5*S/sine
    end do ! end i
    curv(l) = curv(l) + 1/kappa !kappa!1/kappa
end do !end l

!print *, cos_lim1, cos_lim2
hard = 1.

do l = 1, n_chain
    !!!
    !!! n_mon+(l-1)*n_mon is the last bead of chain l.
    !!! We are triggering the active force using this bead.
    !!! It is the most "non-arbitrary" way of setting this
    !!! condition. We could use other beads, but they would 
    !!! excite higher modes of the chain. In this work we 
    !!! need the lowest frequency possible, so we picked 
    !!! the last one. Further study of this modes excitation
    !!! is needed.
    !!! 
    r_2rel = r0(:, n_mon+(l-1)*n_mon)- r0(:, 1+(l-1)*n_mon)  
    cos_th = r_2rel(1)/SQRT(DOT_PRODUCT(r_2rel,r_2rel))
    !k_old = k0(1+(l-1)*(n_mon-1))  ! storing old k

    if( cos_th.gt.cos_lim1 ) then  ! the "more horizontal" (1) section
        !print *, 'is in (1)'
        k0(1+(l-1)*(n_mon-1)) = k_bend
        !do i = 2, n_mon-1
        !    k0(i+(l-1)*(n_mon-1)) = k_bend                
        !end do
    else if( cos_th.lt.cos_lim2 ) then ! the "more vertical" (2) section 
        !print *, 'is in (2)'
        k0(1+(l-1)*(n_mon-1)) = k_bend
        !do i = 2, n_mon-1
        !    k0(i+(l-1)*(n_mon-1)) = k_bend                
        !end do
    else ! the middle section
        ! print *, 'middle and'
        if( cos_mem(l).le.cos_lim2 ) then ! if it came from 2

            kick(l) = kick(l) + 15
            k0(1+(l-1)*(n_mon-1)) = k_act !-.005*k_bend  ! soft k

        else if( cos_mem(l).gt.cos_lim1 ) then ! if it came from 1

            k0(1+(l-1)*(n_mon-1)) = k_bend   !

        end if
    end if
    ! write(55,*) cos_th, cos_lim1, cos_lim2
    cos_mem(l) = cos_th  !saving old position
end do

case(4)
    ! This will be a completely active case

end select

count_local = count_local + 1

if(count_local .eq.10) then
    write(79,'(I12.4)', advance='no') i_time
    write(81,'(I12.4)', advance='no') i_time
    do l=1, n_chain
        write(79,'(ES12.4)', advance='no') kick(l)
        write(81,'(ES12.4)', advance='no') curv(l)
    end do
    write(55,*) i_time, cos_th
    write(79,*) 
    write(81,*) 

    count_local = 1
    kick = 0
    curv = 0.
end if

CONTAINS

FUNCTION cross(a, b) 
    REAL (kind=8), INTENT (in) :: a(3), b(3)
    REAL (kind=8), DIMENSION(3) :: cross
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

!FUNCTION norm1(a) 
!    REAL (kind=8), INTENT (in) :: a(:)
!    REAL (kind=8) :: norm, i, suma  
!    norm = sqrt(dot_product(a,a))
!    !suma = 0
!    !do i = 1, SIZE(a)
!    !    suma = suma + a(i)**2
!    !end do
!    !norm = SQRT(suma)
!END FUNCTION norm
!FUNCTION norm333(a) 
!    REAL (kind=8), INTENT (in) :: a(:)
!    REAL (kind=8) :: norm, i, suma  
!    norm = sqrt(dot_product(a,a))
!    !suma = 0
!    !do i = 1, SIZE(a)
!    !    suma = suma + a(i)**2
!    !end do
!    !norm = SQRT(suma)
!END FUNCTION norm
#endif /*close ACTIVE_BRUSH*/

end subroutine metronome

