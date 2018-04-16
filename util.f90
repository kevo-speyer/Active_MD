module util
#include 'control_simulation.h'
!use commons, only: r0,z_space_wall,histo_b,histo_d,     &
!                n_mon,n_chain,n_mon_tot,n_mon_d,n_chain_d,n_part,i_time &
!                ,n_relax,count_obs
        use commons
        implicit none

contains

subroutine write_conf(mode,vector,length)
!write xyz files for viewing the system in different ways
! mode = 1 : uses 4 different colors
! mode !=1 : uses only one color
! vector: positions of the atoms
! length: length (number of monomers )of the chain to change the colors
implicit none
real (kind=8), dimension(:,:)  :: vector
integer,intent(in) :: mode,length 
!character (len=4) ,dimension(size(vector,dim=1))  :: at_names 
integer :: n,i,ic 
character (len=4) , parameter , dimension(4) :: nam1= (/"Cl  ","N   ","S   ","O   "/)
n=size(vector,dim=1)

open(unit=25,file="config.xyz",status="unknown")


select case (mode)
case(1)

        write(25,*) n
        write(25,*) 'mode=1, long system of 4 colors. chain length = ',length
        ic =1
        do i = 1,n
        if(mod(i-1,length).eq.0) ic = ic + 1
        if(ic > size(nam1) ) ic = 1
        write(25,'(a4,3f10.4)') nam1(ic),vector(i,:)
        end do 
case default
! default uses a unique long conf of similar atoms
        write(25,*) n
        write(25,*) 'mode= default, long system of same color'
        write(25,'("Cl   ",3f10.4)') (vector(i,:),i=1,n)

end select
!write
        close(unit=25)

end subroutine write_conf
!-------------------------------------------------------------------

!    Histograms, for exmaple for monomer profiles
!   mode = 1 : initialize counting
!   mode !=1  : makes histogram in a cummulative way
!   n_step = number of columns in the histogram
!   r_interval = interval length
!   vector: 

!NOW IN functions.f90 subroutine histo(mode) ! ,vector,r_interval,histo_clau)
!NOW IN functions.f90 
!NOW IN functions.f90 !use commons
!NOW IN functions.f90 implicit none
!NOW IN functions.f90 integer,intent(in) :: mode
!NOW IN functions.f90 !real (kind=8), dimension(:), intent(in) :: vector
!NOW IN functions.f90 !real (kind=8), intent(in) :: r_interval
!NOW IN functions.f90 !real (kind=8), dimension(:), intent(inout) :: histo_clau
!NOW IN functions.f90 real (kind=8) :: binbox_c,binstep_c
!NOW IN functions.f90 integer :: n_max,i_step,i_part,n_step_clau
!NOW IN functions.f90 real (kind=8) , pointer ::  vector(:)
!NOW IN functions.f90 real (kind=8) :: r_interval
!NOW IN functions.f90 !real (kind=8) , dimension(size(histo_b) :: histo_clau
!NOW IN functions.f90 
!NOW IN functions.f90 !
!NOW IN functions.f90 
!NOW IN functions.f90  select case (mode)
!NOW IN functions.f90 ! case(0) ! initialization
!NOW IN functions.f90      
!NOW IN functions.f90  case(1) !histogram for polymer brush    
!NOW IN functions.f90 
!NOW IN functions.f90 !     vector => r0(3,1:n_mon*n_chain)
!NOW IN functions.f90 !     histo_clau => histo_b(:)
!NOW IN functions.f90 !    n_step_clau=size(histo_clau,dim=1) 
!NOW IN functions.f90 !    n_max=size(vector,dim=1)  !numebr of particles
!NOW IN functions.f90      r_interval = z_space_wall
!NOW IN functions.f90      vector =>  r0(3,1:n_mon*n_chain)
!NOW IN functions.f90      n_step_clau=size(histo_b,dim=1) 
!NOW IN functions.f90      n_max=size(vector,dim=1)  !numebr of particles
!NOW IN functions.f90 
!NOW IN functions.f90        binstep_c=r_interval/real(n_step_clau)
!NOW IN functions.f90        binbox_c=0.0
!NOW IN functions.f90        do i_step=1,n_step_clau
!NOW IN functions.f90             binbox_c=binbox_c+binstep_c
!NOW IN functions.f90             do i_part=1,n_max
!NOW IN functions.f90                if(vector(i_part).lt.binbox_c) then
!NOW IN functions.f90                    if(vector(i_part).ge.(binbox_c-binstep_c)) then
!NOW IN functions.f90                    histo_b(i_step)=histo_b(i_step)+1
!NOW IN functions.f90                    end if
!NOW IN functions.f90                end if
!NOW IN functions.f90             end do
!NOW IN functions.f90        end do
!NOW IN functions.f90        
!NOW IN functions.f90  case(2) ! histogram for droplet/melt
!NOW IN functions.f90 
!NOW IN functions.f90       r_interval= z_space_wall
!NOW IN functions.f90       vector =>  r0(3,n_mon*n_chain+1:n_mon_tot)
!NOW IN functions.f90       n_step_clau=size(histo_d,dim=1) 
!NOW IN functions.f90       n_max=size(vector,dim=1)  !numebr of particles
!NOW IN functions.f90       binstep_c=r_interval/dble(n_step_clau)
!NOW IN functions.f90       binbox_c=0.0
!NOW IN functions.f90 !      
!NOW IN functions.f90        do i_step=1,n_step_clau
!NOW IN functions.f90             binbox_c=binbox_c+binstep_c
!NOW IN functions.f90             do i_part=1,n_max
!NOW IN functions.f90                if(vector(i_part).lt.binbox_c) then
!NOW IN functions.f90                    if(vector(i_part).ge.(binbox_c-binstep_c)) then
!NOW IN functions.f90                    histo_d(i_step)=histo_d(i_step)+1
!NOW IN functions.f90                    end if
!NOW IN functions.f90                end if
!NOW IN functions.f90             end do
!NOW IN functions.f90        end do
!NOW IN functions.f90 
!NOW IN functions.f90 !       print*,"melt=",sum(histo_d(:),dim=1)/(dble(i_time-n_relax))
!NOW IN functions.f90 end select
!NOW IN functions.f90 
!NOW IN functions.f90 
!NOW IN functions.f90             
!NOW IN functions.f90 
!NOW IN functions.f90 end subroutine histo
! ---------------------------------------------------------
! WARNING: I keep this here but is not being used now.
! for gaussian number generator the module ziggurat 
! is currently in use
!
! ******* Gaussian Random number generator ******
!
! This version is made for using with the r250
! packages. Initialization of them must also be made
! from the calling program.


!     ALGORITHM 488 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12,
!     P. 704.                                                           
!!!REAL FUNCTION GRAND(N)                                                 
!!!! EXCEPT ON THE FIRST CALL GRAND RETURNS A
!!!! PSEUDO-RANDOM NUMBER HAVING A GAUSSIAN (I.E.
!!!! NORMAL) DISTRIBUTION WITH ZERO MEAN AND UNIT
!!!! STANDARD DEVIATION.  THUS, THE DENSITY IS  F(X) =
!!!! EXP(-0.5*X**2)/SQRT(2.0*PI). THE FIRST CALL
!!!! INITIALIZES GRAND AND RETURNS ZERO.
!!!! THE PARAMETER N IS DUMMY.
!!!! GRAND CALLS A FUNCTION RAND, AND IT IS ASSUMED THAT
!!!! SUCCESSIVE CALLS TO RAND(0) GIVE INDEPENDENT
!!!! PSEUDO- RANDOM NUMBERS DISTRIBUTED UNIFORMLY ON (0,
!!!! 1), POSSIBLY INCLUDING 0 (BUT NOT 1).
!!!! THE METHOD USED WAS SUGGESTED BY VON NEUMANN, AND
!!!! IMPROVED BY FORSYTHE, AHRENS, DIETER AND BRENT.
!!!! ON THE AVERAGE THERE ARE 1.37746 CALLS OF RAND FOR
!!!! EACH CALL OF GRAND.
!!!! WARNING - DIMENSION AND DATA STATEMENTS BELOW ARE
!!!!           MACHINE-DEPENDENT.
!!!! DIMENSION OF D MUST BE AT LEAST THE NUMBER OF BITS
!!!! IN THE FRACTION OF A FLOATING-POINT NUMBER.
!!!! THUS, ON MOST MACHINES THE DATA STATEMENT BELOW
!!!! CAN BE TRUNCATED.
!!!! IF THE INTEGRAL OF SQRT(2.0/PI)*EXP(-0.5*X**2) FROM
!!!! A(I) TO INFINITY IS 2**(-I), THEN D(I) = A(I) -
!!!! A(I-1).  
!!!      use commons, only: mz,kptr
!!!      real :: A,W,V,rand(2)
!!!      real :: U,D
!!!      DIMENSION D(60)
!!!      integer :: I,N
!!!      DATA D(1), D(2), D(3), D(4), D(5), D(6), D(7),                    &
!!!     & D(8), D(9), D(10), D(11), D(12), D(13),                          &
!!!     & D(14), D(15), D(16), D(17), D(18), D(19),                        &
!!!     & D(20), D(21), D(22), D(23), D(24), D(25),                        &
!!!     & D(26), D(27), D(28), D(29), D(30), D(31),                        &
!!!     & D(32) /0.674489750,0.475859630,0.383771164,                      &
!!!     & 0.328611323,0.291142827,0.263684322,                             &
!!!     & 0.242508452,0.225567444,0.211634166,                             &
!!!     & 0.199924267,0.189910758,0.181225181,                             &
!!!     & 0.173601400,0.166841909,0.160796729,                             &
!!!     & 0.155349717,0.150409384,0.145902577,                             &
!!!     & 0.141770033,0.137963174,0.134441762,                             &
!!!     & 0.131172150,0.128125965,0.125279090,                             &
!!!     & 0.122610883,0.120103560,0.117741707,                             &
!!!     & 0.115511892,0.113402349,0.111402720,                             &
!!!     & 0.109503852,0.107697617/
!!!      DATA D(33), D(34), D(35), D(36), D(37), D(38),                    &
!!!     & D(39), D(40), D(41), D(42), D(43), D(44),                        &
!!!     & D(45), D(46), D(47), D(48), D(49), D(50),                        &
!!!     & D(51), D(52), D(53), D(54), D(55), D(56),                        &
!!!     & D(57), D(58), D(59), D(60)                                       &
!!!     & /0.105976772,0.104334841,0.102766012,                            &
!!!     & 0.101265052,0.099827234,0.098448282,                             &
!!!     & 0.097124309,0.095851778,0.094627461,                             &
!!!     & 0.093448407,0.092311909,0.091215482,                             &
!!!     & 0.090156838,0.089133867,0.088144619,                             &
!!!     & 0.087187293,0.086260215,0.085361834,                             &
!!!     & 0.084490706,0.083645487,0.082824924,                             &
!!!     & 0.082027847,0.081253162,0.080499844,                             &
!!!     & 0.079766932,0.079053527,0.078358781,                             &
!!!     & 0.077681899/
!!!! END OF MACHINE-DEPENDENT STATEMENTS
!!!! U MUST BE PRESERVED BETWEEN CALLS.
!!!      DATA U /0.0/
!!!      
!!!! INITIALIZE DISPLACEMENT A AND COUNTER I.
!!!      A = 0.0
!!!      I = 0
!!!! INCREMENT COUNTER AND DISPLACEMENT IF LEADING BIT
!!!! OF U IS ONE.
!!!   10 U = U + U
!!!      IF (U.LT.1.0) GO TO 20
!!!      U = U - 1.0
!!!      I = I + 1
!!!      A = A - D(I)
!!!      GO TO 10
!!!! FORM W UNIFORM ON 0 .LE. W .LT. D(I+1) FROM U.
!!!   20 W = D(I+1)*U
!!!! FORM V = 0.5*((W-A)**2 - A**2). NOTE THAT 0 .LE. V
!!!! .LT. LOG(2).
!!!      V = W*(0.5*W-A)
!!!! GENERATE NEW UNIFORM U.
!!!        call r250(mz,rand,2,2,kptr)
!!!   30 U = rand(1) ! ori RAND(0)
!!!! ACCEPT W AS A RANDOM SAMPLE IF V .LE. U.
!!!      IF (V.LE.U) GO TO 40
!!!! GENERATE RANDOM V.
!!!      V = rand(2) ! ori RAND(0)
!!!! LOOP IF U .GT. V.
!!!      IF (U.GT.V) GO TO 30
!!!! REJECT W AND FORM A NEW UNIFORM U FROM V AND U.
!!!      U = (V-U)/(1.0-U)
!!!      GO TO 20
!!!! FORM NEW U (TO BE USED ON NEXT CALL) FROM U AND V.
!!!   40 U = (U-V)/(1.0-V)
!!!! USE FIRST BIT OF U FOR SIGN, RETURN NORMAL VARIATE.
!!!      U = U + U
!!!      IF (U.LT.1.0) GO TO 50
!!!      U = U - 1.0
!!!      GRAND = W - A
!!!      RETURN
!!!   50 GRAND = A - W
!!!      RETURN
!!!        END FUNCTION GRAND 
! --
! ********* General histogram generation with a vector as input *********
subroutine histo_vec(vector) ! ,vector,r_interval,histo_clau)
!use commons
implicit none
real , dimension(:), intent(in) :: vector
!real (kind=8), dimension(:), intent(in) :: vector
real (kind=8)  :: r_interval
real (kind=8) :: binbox_c,binstep_c
integer :: n_max,i_step,i_part,n_step_clau
real (kind=8) :: histo_v(200),eje_x(200)
!real (kind=8) , dimension(size(histo_v) :: histo_clau
r_interval = maxval(vector) - minval(vector)
n_step_clau=size(histo_v,dim=1)
n_max=size(vector,dim=1)  !numebr of particles

binstep_c=r_interval/real(n_step_clau)
binbox_c= minval(vector)
do i_step=1,n_step_clau
binbox_c=binbox_c+binstep_c
eje_x(i_step) = binbox_c
do i_part=1,n_max
if(vector(i_part).lt.binbox_c) then
if(vector(i_part).ge.(binbox_c-binstep_c)) then
histo_v(i_step)=histo_v(i_step)+1
end if
end if
end do
end do

do i_step = 1 , n_step_clau
print *,eje_x(i_step),histo_v(i_step)/real(n_max)
end do
end subroutine histo_vec
!--------------------------------------

      subroutine check_fluid_fluid
      use commons
      implicit none


      
      v_fluid_fluid = 0.

!
 do i_part = 1,n_part-1
       i_type =a_type(i_part)
       do j_part = i_part+1, n_part
        j_type =a_type(j_part)

         delta_r(1) = r0(1,i_part) - r0(1,j_part)
         delta_r(2) = r0(2,i_part) - r0(2,j_part)
         delta_r(3) = r0(3,i_part) - r0(3,j_part)

!-----  get boundaries right
#   if SYMMETRY == 0 
        do i_dim = 1,n_dim  -1
#   elif SYMMETRY == 1            
        do i_dim = 1,n_dim 
#endif

         delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)*inv_boundary(i_dim))
        end do
        r_2 = 0.
        do i_dim = 1,n_dim
         r_2 = r_2 + delta_r(i_dim)**2
        end do
!-----  check whether interaction takes place
      if(r_2.lt.range_2(i_type,j_type)) then
         r_2 = max(r_2,r_2_min)!(i_type, j_type))
         r_6 = (sigma_2(i_type,j_type)/r_2)**3
         r_12 = r_6**2
         pot_loc = (r_12-r_6) - e_shift(i_type,j_type)
         v_fluid_fluid = v_fluid_fluid + epsil(i_type,j_type)*pot_loc
         r_dummy = epsil(i_type,j_type)*(-12*r_12+6*r_6)/r_2
! NOTE: no force calculation
      end if
!
      end do ! loop over particles
end do ! loop over particles
 print *,"[check_fluid_fluid] V= ",v_fluid_fluid
          !
end subroutine check_fluid_fluid

! ------------------------------------
! Velocity profile 

! Calculates a mean velocity profile along z coordinate
! mode : routine mode
! z: positions, only z coordinate
! v: velocities
! z_width: sample width
! n_layers: number of layers in which the profile will be done
! velocity profile

!NOW IN velocity_prof subroutine velocity_prof(mode,z_min,z_width)
!NOW IN velocity_prof implicit none
!NOW IN velocity_prof integer ,intent(in) :: mode
!NOW IN velocity_prof !real (kind=8) , dimension(:,:), intent(in) :: v
!NOW IN velocity_prof !real (kind=8) , dimension(:)  , intent(in) :: z
!NOW IN velocity_prof real (kind=8) , intent(in) , optional   :: z_width,z_min
!NOW IN velocity_prof !real (kind=8) , dimension(n_layers,3), save :: v_prof_b,v_prof_b_2,v_prof_m,v_prof_m_2
!NOW IN velocity_prof real (kind=8) , dimension(n_layers,3):: v_prof_b,v_prof_b_2,v_prof_m,v_prof_m_2
!NOW IN velocity_prof 
!NOW IN velocity_prof integer :: i,j,ii,jj,kk,z_count_b,z_count_m
!NOW IN velocity_prof integer, save :: n_parti
!NOW IN velocity_prof integer       :: n_time = 0
!NOW IN velocity_prof real(kind=8) :: zmin,zmax,v_cell_b(3),v_cell_m(3)
!NOW IN velocity_prof real(kind=8) ,save :: dz 
!NOW IN velocity_prof 
!NOW IN velocity_prof select case (mode)
!NOW IN velocity_prof 
!NOW IN velocity_prof         case(1) ! Initialize 
!NOW IN velocity_prof             v_prof_b(:,:) = 0.
!NOW IN velocity_prof             v_prof_b_2(:,:) = 0.
!NOW IN velocity_prof             v_prof_m(:,:) = 0.
!NOW IN velocity_prof             v_prof_m_2(:,:) = 0.
!NOW IN velocity_prof             n_time = 0 
!NOW IN velocity_prof             dz = z_width / real(n_layers)
!NOW IN velocity_prof             n_parti = size(r0(3,:))
!NOW IN velocity_prof             print*,"*  Velocity prof initialized"
!NOW IN velocity_prof         case(2) ! profile calculation
!NOW IN velocity_prof !            print *,n_time,dz
!NOW IN velocity_prof             n_time = n_time + 1
!NOW IN velocity_prof             do i = 1,n_layers
!NOW IN velocity_prof             v_cell_b(:) = 0.
!NOW IN velocity_prof             v_cell_m(:) = 0.
!NOW IN velocity_prof             z_count_b = 0
!NOW IN velocity_prof             z_count_m = 0
!NOW IN velocity_prof             zmin = z_min + dz*(real(i) - 1.)
!NOW IN velocity_prof             zmax = z_min + dz*real(i)
!NOW IN velocity_prof                do j= 1,n_parti 
!NOW IN velocity_prof                  if( r0(3,j) < zmax.and.r0(3,j) > zmin) then
!NOW IN velocity_prof                      
!NOW IN velocity_prof                    if(a_type(j) == 2 ) then  ! If particle belongs to the brush
!NOW IN velocity_prof                    z_count_b = z_count_b + 1    
!NOW IN velocity_prof                    v_cell_b(:) = v_cell_b(:) + dt*v(:,j)/dt
!NOW IN velocity_prof                    end if
!NOW IN velocity_prof                    
!NOW IN velocity_prof                    if(a_type(j) == 3 ) then  ! If particle belongs to the melt 
!NOW IN velocity_prof                    z_count_m = z_count_m + 1    
!NOW IN velocity_prof                    v_cell_m(:) = v_cell_m(:) + dt*v(:,j)/dt
!NOW IN velocity_prof                    end if
!NOW IN velocity_prof 
!NOW IN velocity_prof                  end if
!NOW IN velocity_prof                end do 
!NOW IN velocity_prof ! WARNING: probe not renormalizing the velocities. (I can then divide de density profile)                 
!NOW IN velocity_prof                  if(z_count_b > 0 ) then
!NOW IN velocity_prof                  v_prof_b(i,:) = v_prof_b(i,:) + v_cell_b(:)                ! / real(z_count_b)
!NOW IN velocity_prof                  v_prof_b_2(i,:) = v_prof_b_2(i,:) + (v_cell_b(:))**2       !/ real(z_count_b)
!NOW IN velocity_prof                  end if
!NOW IN velocity_prof 
!NOW IN velocity_prof                  if(z_count_m > 0 ) then
!NOW IN velocity_prof                  v_prof_m(i,:) = v_prof_m(i,:) + v_cell_m(:)                !/ real(z_count_m)
!NOW IN velocity_prof                  v_prof_m_2(i,:) = v_prof_m_2(i,:) + (v_cell_m(:))**2       !/ real(z_count_m)
!NOW IN velocity_prof                  end if
!NOW IN velocity_prof             end do
!NOW IN velocity_prof !   Debug            
!NOW IN velocity_prof !            do i = 1,n_parti
!NOW IN velocity_prof !            print '(12f16.3)',v_prof_b(1,:),v_prof_b_2(1,:)
!NOW IN velocity_prof !            end do 
!NOW IN velocity_prof         case(3) ! Final normalization 
!NOW IN velocity_prof             v_prof_b(:,:) = v_prof_b(:,:) / real(n_time)
!NOW IN velocity_prof             v_prof_b_2(:,:) = v_prof_b_2(:,:) / real(n_time)
!NOW IN velocity_prof             v_prof_m(:,:) = v_prof_m(:,:) / real(n_time)
!NOW IN velocity_prof             v_prof_m_2(:,:) = v_prof_m_2(:,:) / real(n_time)
!NOW IN velocity_prof             open(unit=74)
!NOW IN velocity_prof             write(74,'(13a15)') "#z","vx_brush","vy_brush","vz_brush","vx2_brush","vy2_brush","vz2_brush" &
!NOW IN velocity_prof                                 ,"vx_melt","vy_melt","vz_melt","vx2_melt","vy2_melt","vz2_melt"
!NOW IN velocity_prof             do i = 1 , n_layers
!NOW IN velocity_prof             write(74,'(13f17.6)') z_min+dz/2.+real(i-1)*dz,v_prof_b(i,:),v_prof_b_2(i,:),                       &
!NOW IN velocity_prof                                                      v_prof_m(i,:),v_prof_m_2(i,:)
!NOW IN velocity_prof             end do
!NOW IN velocity_prof         case default
!NOW IN velocity_prof             print*,"*mode* must have the values 1,2, or 3. Stop here."
!NOW IN velocity_prof             stop
!NOW IN velocity_prof 
!NOW IN velocity_prof end select            
!NOW IN velocity_prof end subroutine velocity_prof


end module util

