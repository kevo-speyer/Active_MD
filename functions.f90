module functions
        use commons
#include "control_simulation.h"                
        implicit none
!        integer :: i,j,k,ii,jj,kk ! general counters variables

contains
! Use really these routines if PROFILES is on.

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

!  **** Periodic Boudary conditions refolding
!
        subroutine refold(mode,n_mon_loc,n_ini,n_end)
!        real(kind=8) ,intent(in) :: 
!        real(kind=8) ,intent(out) :: spabs(:,:) 
        integer , intent(in) :: mode,n_mon_loc,n_ini, n_end
        integer , save :: n_part , n_dim_loc
        integer :: i_dim,i_part
        select case(mode)
        case(0)
#   if SYMMETRY == 0         
                n_dim_loc = 2 
#   elif SYMMETRY == 1 /* bulki */
                n_dim_loc = 3 
#   endif 
        case(1)
              r0_unfold(:,n_ini:n_end) = r0(:,n_ini:n_end)
              do i_part=n_ini,n_end
                 do i_dim=1,n_dim_loc
                   if(mod(i_part-1,n_mon_loc).ne.0) then
                      if(r0_unfold(i_dim,i_part)-r0_unfold(i_dim,i_part-1).gt.half_boundary(i_dim)) then
                          r0_unfold(i_dim,i_part)=r0_unfold(i_dim,i_part)-boundary(i_dim)
                          cycle
                      end if
                      if(r0_unfold(i_dim,i_part-1)-r0_unfold(i_dim,i_part).gt.half_boundary(i_dim)) then
                          r0_unfold(i_dim,i_part)=r0_unfold(i_dim,i_part)+boundary(i_dim)
                          cycle
                      end if
                   end if
                 end do
              end do
          end select     

          end subroutine refold

#ifdef PROFILES
!-------------------------------------------------------------------

!    Histograms, for monomer profiles
!   mode = 0 : initialize counting
!   mode = 1  : makes histogram  for brushes
!   mode = 2  : makes histogram  for melt
!   mode = 3  : writes out final data
!   
subroutine histo(mode,hist_dim)

    implicit none
    real (kind=8) ,allocatable ,save:: dens_z(:,:),r_scaled(:,:)
    real (kind=8) :: r_bin
    real (kind=8), save :: z_step,inv_z_step,surf
    integer,intent(in) :: mode,hist_dim
    integer ,save :: count_obs_b,count_obs_m
    integer :: n_max,i_step,i_part,n_step_clau,i_dummy
    integer :: z_index

 select case (mode)
 case(0) ! initialization

!  histo(hist_dim,1)  : top brush || brush, if notdef BRUSH_SEP
!  histo(hist_dim,2)  : bottom brush
!  histo(hist_dim,2) : melt 
!  histo(hist_dim,3) : particle 4 

#ifdef BRUSH_SEP
    allocate(dens_z(hist_dim,4),r_scaled(3,n_mon_tot))

#else

    allocate(dens_z(hist_dim,3),r_scaled(3,n_mon_tot))
#endif

    r_scaled = 0.

    dens_z(:,:) = 0.
    z_step = boundary(3) / dble(hist_dim)  
    inv_z_step = 1./z_step

    surf = boundary(1)*boundary(2)

    count_obs_b = 0
!ori    count_obs_m = 0

     case(1)  ! Writing files and normalizing

        print '(/a,i6/)','   *  Density profiles over steps= ',count_obs_b 

        open(unit=73,file='dens_prof.mide',status='unknown')


         dens_z(:,:) =  dens_z(:,:) * dble(hist_dim)*inv_boundary(3)/ (dble(count_obs_b)*surf)

         histo_out(:,:) = dens_z(:,:)  ! global: norm_vel prof and sigma prof are calculated with this vector 

         r_dummy = boundary(3)/dble(hist_dim)
         r_bin = r_dummy/2.
         do i_dummy = 1,hist_dim

#ifdef BRUSH_SEP
            write(73,'(5f16.7)') r_bin,dens_z(i_dummy,1),dens_z(i_dummy,2),dens_z(i_dummy,3), dens_z(i_dummy,4)
#else

            write(73,'(4f16.7)') r_bin,dens_z(i_dummy,1),dens_z(i_dummy,2),dens_z(i_dummy,3)
#endif
            r_bin = r_bin + r_dummy
         end do
          

        case(2)   ! Order N algorithm

        count_obs_b = count_obs_b + 1

         r_scaled(3,:) = r0(3,:)*inv_z_step

!        do i_part = 1,n_mon_tot
!        if(r0(3,i_part) > 30. ) print *i_time 
!        end do
! Brush

#ifdef BRUSH_SEP
        do i_part = 1, part_init_d/2
#else
        do i_part = 1, part_init_d
#endif
            z_index = int (r_scaled(3,i_part) ) + 1

            if (z_index<1) z_index = 1
            if (z_index>hist_dim) z_index = hist_dim

            dens_z(z_index,1) =  dens_z(z_index,1) + 1.
        end do
! Melt / Droplet
        do i_part = part_init_d+1,part_init_e
            z_index = int (r_scaled(3,i_part) ) + 1

            if (z_index<1) z_index = 1
            if (z_index>hist_dim) z_index = hist_dim

            dens_z(z_index,2) =  dens_z(z_index,2) + 1.
        end do

! Particle 4 if applies
#ifdef PARTICLE_4
        do i_part = part_init_e+1,n_mon_tot

            z_index = int (r_scaled(3,i_part) ) + 1

            if (z_index<1) z_index = 1
            if (z_index>hist_dim) z_index = hist_dim

            dens_z(z_index,3) =  dens_z(z_index,3) + 1.

        end do
#endif
#ifdef BRUSH_SEP
do i_part = part_init_d/2+1,part_init_d
            z_index = int (r_scaled(3,i_part) ) + 1

            if (z_index<1) z_index = 1
            if (z_index>hist_dim) z_index = hist_dim

            dens_z(z_index,4) =  dens_z(z_index,4) + 1.
        end do

#endif
end select
end subroutine histo
!
! ******** histo_vec *********
!
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
!
! **** Velocity profile 

! Calculates a mean velocity profile along z coordinate
! mode : routine mode
! z: positions, only z coordinate
! v: velocities
! z_width: sample width
! n_layers: number of layers in which the profile will be done
! velocity profile

#ifdef PROFILES
subroutine velocity_prof(mode,n_layers,z_min,z_width) ! ,v_prof,v_prof_2) ! NOW in commons 
    implicit none
    integer ,intent(in) :: mode,n_layers
    real (kind=8) , intent(in) , optional   :: z_width,z_min
    integer :: i,j,ii,jj,kk,z_count_b,z_count_m
    real(kind=8) :: zmin,zmax,v_cell_b(3),v_cell_m(3)
    real(kind=8) , allocatable, save ::  r_scaled(:,:)
    real(kind=8) ,save :: dz,inv_z_step
    real(kind=8)       :: inv_n_time
    integer, save :: n_time
    integer       ::  z_index

select case (mode)

    case(0) ! Initialize 

    v_prof(:,:,:) = 0.
    v_prof_2(:,:,:) = 0.

    allocate( r_scaled(3,n_mon_tot) )

    !            n_time = 0 
    dz = boundary(3) /dble(n_layers)
    inv_z_step = 1/dz

    print '(/a/)',"   *  Velocity prof initialized. "

    n_time = 0

    case(1) ! Final normalization 

     
    inv_n_time = 1/dble(n_time)

    v_prof(:,:,:)   = v_prof(:,:,:)   *inv_n_time 
    v_prof_2(:,:,:) = v_prof_2(:,:,:) *inv_n_time    


! velocity profile

    open(unit=74,file='vel_prof.mide',status='unknown')
#ifndef PARTICLE_4
    write(74,'(13a15)') "#z","vx_brush","vy_brush","vz_brush","vx2_brush","vy2_brush","vz2_brush" &
    ,"vx_melt","vy_melt","vz_melt","vx2_melt","vy2_melt","vz2_melt"
    do i = 1 , n_layers
    write(74,'(13e16.7)') z_min+dz/2.+dble(i-1)*dz,v_prof(i,:,1),v_prof_2(i,:,1),v_prof(i,:,2),v_prof_2(i,:,2)
    end do
#else
    write(74,'(19a16)') "#z","vx_brush","vy_brush","vz_brush","vx2_brush","vy2_brush","vz2_brush" &
    ,"vx_melt","vy_melt","vz_melt","vx2_melt","vy2_melt","vz2_melt", "vx_part4","vy_part4","vz_part4","vx2_part4","vy2_part4","vz2_part4"
    do i = 1 , n_layers
    write(74,'(19e16.7)') z_min+dz/2.+dble(i-1)*dz,v_prof(i,:,1),v_prof_2(i,:,1),v_prof(i,:,2),v_prof_2(i,:,2), &
                                                   v_prof(i,:,3),v_prof_2(i,:,3)
    end do
#endif

    case(2) ! Profile calculation

         n_time = n_time + 1

         r_scaled(3,:) = r0(3,:)*inv_z_step

! Brush
        do i_part = 1, part_init_d

            z_index = int (r_scaled(3,i_part) ) + 1

            if(z_index<1) z_index =1  
            if(z_index>hist_dim) z_index = hist_dim  

            v_prof(z_index,:,1) =   v_prof(z_index,:,1) + v(:,i_part) 
            v_prof_2(z_index,:,1) =   v_prof_2(z_index,:,1) + v(:,i_part)*v(:,i_part)
        end do

! Melt / Droplet

        do i_part = part_init_d+1,part_init_e

            z_index = int (r_scaled(3,i_part) ) + 1

            if(z_index<1) z_index =1  
            if(z_index>hist_dim) z_index = hist_dim  

            v_prof(z_index,:,2) =   v_prof(z_index,:,2) + v(:,i_part) 
            v_prof_2(z_index,:,2) =   v_prof_2(z_index,:,2) + v(:,i_part)*v(:,i_part)

!        print *,v_prof(:,1,2) ; stop
        
        end do

! Particle 4 if applies
#ifdef PARTICLE_4
        do i_part = part_init_e+1,n_mon_tot

            z_index = int (r_scaled(3,i_part) ) + 1

            if(z_index<1) z_index =1  
            if(z_index>hist_dim) z_index = hist_dim  

            v_prof(z_index,:,3) =   v_prof(z_index,:,3) + v(:,i_part) 
            v_prof_2(z_index,:,3) =   v_prof_2(z_index,:,3) + v(:,i_part)*v(:,i_part)
        end do
#endif

! obsolete: order N^2    !            print *,n_time,dz
! obsolete: order N^2    n_time = n_time + 1
! obsolete: order N^2    do i = 1,n_layers
! obsolete: order N^2    v_cell_b(:) = 0.
! obsolete: order N^2    v_cell_m(:) = 0.
! obsolete: order N^2    z_count_b = 0
! obsolete: order N^2    z_count_m = 0
! obsolete: order N^2    zmin = z_min + dz*(real(i) - 1.)
! obsolete: order N^2    zmax = z_min + dz*real(i)
! obsolete: order N^2    do j= 1,n_parti 
! obsolete: order N^2
! obsolete: order N^2    if( r0(3,j) < zmax.and. r0(3,j) > zmin) then
! obsolete: order N^2
! obsolete: order N^2    if(a_type(j) == 2 ) then  ! If particle belongs to the brush
! obsolete: order N^2    z_count_b = z_count_b + 1    
! obsolete: order N^2        v_cell_b(:) = v_cell_b(:) + v(j,:)
! obsolete: order N^2    end if
! obsolete: order N^2
! obsolete: order N^2    if(a_type(j) == 3 ) then  ! If particle belongs to the melt 
! obsolete: order N^2    z_count_m = z_count_m + 1    
! obsolete: order N^2        v_cell_m(:) = v_cell_m(:) + v(j,:)
! obsolete: order N^2    end if
! obsolete: order N^2
! obsolete: order N^2    end if
! obsolete: order N^2    end do 
! obsolete: order N^2
! obsolete: order N^2    v_prof(i,:,1) = v_prof(i,:,1) + v_cell_b(:)   ! / real(z_count_b)
! obsolete: order N^2v_prof_2(i,:,1) = v_prof_2(i,:,1) + (v_cell_b(:))**2  !/ real(z_count_b)
! obsolete: order N^2
! obsolete: order N^2    v_prof(i,:,2) = v_prof(i,:,2) + v_cell_m(:)    ! / real(z_count_m)
! obsolete: order N^2v_prof_2(i,:,2) = v_prof_2(i,:,2) + (v_cell_m(:))**2    ! / real(z_count_m)
! obsolete: order N^2
! obsolete: order N^2    end do  

    case default
    print*,"*mode* must have the values 1,2, or 3. Stop here."
    stop

    end select            
end subroutine velocity_prof
#endif

!
! **** Diffusion coefficient 
! 
#   ifdef DIFF 
    subroutine diff_coef(mode,time)
        integer, intent(in) :: mode
        real (kind=8), intent(in) :: time
        real (kind=8)    :: vec(3),Dm(3),Dt(3)
        integer :: i 
        select case(mode)
        case(1) ! init 


            r0_ref(:,:) = r0_unfold(:,:) 

          print '(/a/)',"   * Initizializing diffusion calculation [by now, only particle 4 ]"
! [Prepared, but not tested]            open(unit=400,file='diff_m.mide',status='unknown')
            open(unit=401,file='diff_t.mide',status='unknown')

 ! [Prepared, but not tested]           write(400,*) "#real_time", " Diffusion coefficient: dx^2 dy^2 dz^2 "
            write(401,*) "#real_time", " Diffusion coefficient: dx^2 dy^2 dz^2 "

        case(2) ! for melt [WARN: has not been tested ] 
! Diffusion coeficient for the melt only             
            Dm(:) = 0. ! D is dr^2
            do i=part_init_d+ 1,part_init_e ! start in the first atom of the melt
            vec(:)= r0_unfold(:,i) - r0_ref(:,i)
            Dm(:) = Dm(:) + vec(:)*vec(:)
            end do
            write(400,'(2f16.7)') time, Dm(:) 

        case(3) ! ---------------- for tagged particles 

             Dt(:) = 0. ! D is dr^2
            do i = part_init_e+ 1,n_mon_tot ! start in the first atom of the melt
            vec(:)= r0_unfold(:,i) - r0_ref(:,i)
            Dt(1) = Dt(1) + vec(1)*vec(1)
            Dt(2) = Dt(2) + vec(2)*vec(2)
            Dt(3) = Dt(3) + vec(3)*vec(3)

            end do
            write(401,'(4f16.7)') time, Dt(:)/dble(n_tot_e) 
        end select
        
        end subroutine diff_coef
#endif /* DIFF */
!        
!                                **** Radius of Gyration  ****
! WARN: assumes equal masses  
        subroutine r_gyr(mode,i_case,n_mon,spabs,r_g2)
        integer , intent(in) :: mode,n_mon,i_case
        real (kind=8) ,intent(in) :: spabs(:,:)
!        real (kind=8) ,intent(out) :: r_g_vec(:,:)
        real (kind=8)  :: r_g2(:,:)
        real (kind=8) , allocatable,save :: r_g2_mean(:,:) 
        integer, allocatable, save :: n_count(:)
        real (kind=8)  :: v_dummy(3) , r_cm(3), r_dummy(3)
        integer :: i_chain, n_chain,i_part,i
        integer, save :: n_case

        select case(mode)
        case(1)  !  init
            n_case = size(r_g2,dim=1)
          allocate( r_g2_mean(n_case,3) , n_count(n_case) ) 
          n_count(:) = 0
          r_g2_mean(:,:) = 0.
        case(2)  !  calculation
              n_count(i_case) = n_count(i_case) + 1
              n_chain = size (spabs,dim = 1) /n_mon 
!              print *,"i_case=",i_case,'n_chain=',n_chain
              r_g2(i_case,:) = 0.d0
!
!        **** Rg^2 **** 
!
                do i_chain = 1,n_chain
                    v_dummy(:) = 0.0
                    r_cm(:) = sum(spabs( (i_chain-1)*n_mon+1:i_chain*n_mon,:),dim=1) / dble(n_mon)

                    do i_part=n_mon*(i_chain-1)+1 ,i_chain*n_mon
                        v_dummy(:) = v_dummy(:) +   (spabs(:,i_part) - r_cm(:) )**2
                    end do
                    v_dummy(:) = v_dummy(:) / dble(n_mon) ! Rg^2 of one polymer
                    r_g2(i_case,:) =  r_g2(i_case,:) + v_dummy(:)  ! mean Rg*N
                end do
                r_g2(i_case,:) =  r_g2(i_case,:)/dble(n_chain)! mean Rg
                r_g2_mean(i_case,:) = r_g2_mean(i_case,:) + r_g2(i_case,:)

        case(3)    
              do i = 1,n_case
              r_g2_mean(i,:) = r_g2_mean(i,:) / real(n_count(i))
              end do
              print '(/15x,a/)', "Radius of Gyration:  "
              print '(a30,3f16.7)',"Top brushes, Rg^2= ",r_g2_mean(1,:)
              print '(a30,3f16.7)',"Bottom brushes, Rg^2= ",r_g2_mean(2,:)
              print '(a30,3f16.7)',"Melt,  Rg^2= ",r_g2_mean(3,:)

              !   r_g2(3,1) = dot_product( r_g2(1,:),r_g2(1,:))
              ! Rg for Melt
!!!!              do i_chain = n_chain+1,n_chain + n_chain_d
!!!!              v_dummy(:) = 0.0
!!!!              r_cm(:) = sum(spabs((i_chain-1)*n_mon_d+1:i_chain*n_mon_d,:),dim=1) / dble(n_mon_d)
!!!!              do i_part= (i_chain-1)*n_mon_d+1,i_chain*n_mon_d
!!!!              v_dummy(:) = v_dummy(:) +   spabs(:,i_part)**2 - r_cm(:)**2
!!!!              end do
!!!!              v_dummy(:) = v_dummy(:) / dble(n_mon_d)! Rg of one polymer
!!!!              r_g2(2,:) =  r_g2(2,:) + v_dummy(:)! sumN Rg
!!!!              end do
!!!!              r_g2(2,:) =  r_g2(2,:)/dble(n_chain_d) ! mean Rg^2
         end select             
!              r_g2_mean(1,:) = r_g2_mean(1,:)+ r_g2(1,:) !brushes
!              r_g2_mean(2,:) = r_g2_mean(2,:)+ r_g2(2,:) ! droplet/melt
!              r_g2_mean(3,:) = r_g2_mean(3,:)+ (r_g2(1,:) + r_g2(2,:))/2.0 ! total
!              
!              r_par_per2(1) = r_par_per2(1)+ (r_g2(1,3)- 0.5*(r_g2(1,1)+r_g2(1,2)))/sum(r_g2(1,:))
!              r_par_per2(2) = r_par_per2(2)+ (r_g2(2,3)- 0.5*(r_g2(2,1)+r_g2(2,2)))/sum(r_g2(2,:))
              ! ------------ end gyration stuff
              !print*, sqrt(dot_product(r_g2(1,:),r_g2(1,:))), sqrt(dot_product(r_g2(2,:),r_g2(2,:)))
              ! energies out

        end subroutine r_gyr


!        
!                                **** Radius of Gyration  2nd Version ****
! WARN: assumes equal masses  
        subroutine r_gyr_2v(mode,n_mon,n_chain,n_mon_d,n_chain_d,spabs,r_g2,r_g_vec,r_cm_vec)
        integer , intent(in) :: mode,n_mon,n_mon_d,n_chain,n_chain_d
        real (kind=8) ,intent(in) :: spabs(:,:)
        real (kind=8) ,intent(out) :: r_g_vec(:,:), r_cm_vec(:,:)
        real (kind=8)  :: r_g2(:,:)
        real (kind=8) , allocatable,save :: r_g2_mean(:,:) 
        integer,  save :: n_count
        real (kind=8)  :: v_dummy(3) , r_cm(3), r_dummy(3)
        integer :: i_chain, i_part,i
        integer, save :: n_case

        select case(mode)
        case(1)  !  init
            n_case = size(r_g2,dim=1)
          allocate( r_g2_mean(n_case,3)  ) 
          n_count = 0
          r_g2_mean(:,:) = 0.

        case(2)  !  calculation
              n_count = n_count + 1

!
!        **** Rg^2 **** 

!  CM calculation  
             do i_chain = 1,n_chain                    ! brush
                  r_cm(:) = sum(spabs( (i_chain-1)*n_mon+1:i_chain*n_mon,:),dim=1) / dble(n_mon)
                  r_cm_vec(i_chain,:) = r_cm(:)
             end do
             do i_chain = n_chain+1, n_chain+n_chain_d ! melt
                  r_cm(:) = sum(spabs( (i_chain-1)*n_mon_d+1:i_chain*n_mon_d,:),dim=1) / dble(n_mon_d)
                  r_cm_vec(i_chain,:) = r_cm(:)
             end do
!  R_gyr itself 
              do i_chain = 1,n_chain   ! brushes
                  v_dummy(:) = 0.0
                  do i_part=n_mon*(i_chain-1)+1 ,i_chain*n_mon ! loop over particles of the same polymer
                  v_dummy(:) = v_dummy(:) +   (spabs(:,i_part) - r_cm_vec(i_chain,:) )**2
                  end do
                  v_dummy(:) = v_dummy(:) / dble(n_mon) ! Rg^2 of one polymer
                  r_g_vec(i_chain,:) = v_dummy(:) 
              end do
              do i_chain = n_chain+1,n_chain+n_chain_d  ! melt 
                  v_dummy(:) = 0.0
                  do i_part=n_mon_d*(i_chain-1)+1 ,i_chain*n_mon_d ! loop over particles of the same polymer
                  v_dummy(:) = v_dummy(:) +   (spabs(:,i_part) - r_cm_vec(i_chain,:) )**2
                  end do
                  v_dummy(:) = v_dummy(:) / dble(n_mon) ! Rg^2 of one polymer
                  r_g_vec(i_chain,:) = v_dummy(:) 
              end do
! [system dependent]              
                  r_g2(1,:) = sum (r_g_vec(1:n_chain/2,:),dim=1) /real(n_chain/2) 
                  r_g2(2,:) = sum (r_g_vec(1+n_chain/2:n_chain,:),dim=1) /real(n_chain/2) 
                  r_g2(3,:) = sum (r_g_vec(1+n_chain:n_chain+n_chain_d,:),dim=1) /real(n_chain_d) 
                  do i = 1,n_case 
                  r_g2_mean(i,:) = r_g2_mean(i,:) + r_g2(i,:)
                  end do
        case(3)    
              do i = 1,n_case
              r_g2_mean(i,:) = r_g2_mean(i,:) / real(n_count)
              end do
              print '(/15x,a/)', "Radius of Gyration:                                         Rg^2  "
              print '(a30,4f16.7)',"Top brushes, Rg^2= ",   r_g2_mean(1,:),sum( r_g2_mean(1,:)) 
              print '(a30,4f16.7)',"Bottom brushes, Rg^2= ",r_g2_mean(2,:),sum( r_g2_mean(2,:))
              print '(a30,4f16.7)',"Melt,  Rg^2= ",         r_g2_mean(3,:),sum( r_g2_mean(3,:))

         end select             

        end subroutine r_gyr_2v
        
!
!
!                                **** End to End Radius  **** SECOND VERSION
!
        subroutine r_end_to_end_2v(mode,n_mon,n_chain,n_mon_d,n_chain_d,r0_unfold,r_ee,r_ee_mean,r_ee2,r_ee_mean2,r_ee_vec)
        integer , intent(in) :: mode,n_mon,n_chain,n_mon_d,n_chain_d
        real (kind=8) , intent(in)   :: r0_unfold(:,:)
        real (kind=8) , intent(out)  :: r_ee(:,:),r_ee_mean(:,:),r_ee2(:,:),r_ee_mean2(:,:),r_ee_vec(:,:) 
        integer,  save :: n_count,n_case,n_chain_tot
        real (kind=8)  :: v_dummy(3) , r_dummy(3),dist
        integer :: i_chain, i_part,i,i_case

        select case(mode)
        case(1)  !  init
            n_case = 3 ! size(r_ee,dim=1)
!          allocate( r_ee_mean(n_case,3) ,n_count(n_case) ) 
          n_count = 0
          r_ee_mean(:,:) = 0.
          r_ee_mean2(:,:) = 0.
          open(unit=76,file='R_ee_dist.mide',status='unknown')
          write(76,*)  '#Distance of Ree for brushes vs time'
        case(2)  !  calculation
              n_count = n_count + 1  
              n_chain_tot = n_chain+n_chain_d   !size (spabs,dim = 1) /n_mon 
!              print*,"n_chain=",n_chain
              r_ee(:,:) = 0.d0
              r_ee2(:,:) = 0.d0
!
!        **** R end to end = Ree **** 
!
              do i_chain = 1,n_chain_tot
                 v_dummy(:) = spabs(i_chain*n_mon, :) - spabs((i_chain-1)*n_mon+1,:) 
          ! for each chain       
                  r_ee_vec(i_chain,:) = v_dummy(:)
                  end do
!                  
! dintiguish top,  bottom wall and melt for mean values
! [this might be system dependent ]         

                 do i_chain=1,n_chain/2                     ! top
                  r_ee(1,:)=r_ee(1,:)+ r_ee_vec(i_chain,:)
                  r_ee2(1,:)=r_ee2(1,:)+ r_ee_vec(i_chain,:)**2 
                  dist = sqrt( dot_product (r_ee_vec(i_chain,:),r_ee_vec(i_chain,:) ) ) 
                 write(76,*) n_count,dist
                 end do
                  
                 do i_chain=n_chain/2+1 ,n_chain            ! bottom
                  r_ee(2,:)=r_ee(2,:)+ r_ee_vec(i_chain,:)
                  r_ee2(2,:)=r_ee2(2,:)+ r_ee_vec(i_chain,:)**2 
                 dist = sqrt( dot_product( v_dummy(:),v_dummy(:) ) ) 
                 write(76,*) n_count,dist
                 end do

                 do i_chain=n_chain +1, n_chain+n_chain_d   ! melt
                  r_ee(3,:)=r_ee(3,:)+ r_ee_vec(i_chain,:)
                  r_ee2(3,:)=r_ee2(3,:)+ r_ee_vec(i_chain,:)**2 
                 end do
                 
              r_ee(1,:) =  r_ee(1,:)/dble(n_chain/2)
              r_ee(2,:) =  r_ee(2,:)/dble(n_chain/2)
              r_ee(3,:) =  r_ee(3,:)/dble(n_chain_d)

              r_ee2(1,:) =  r_ee2(1,:)/dble(n_chain/2)
              r_ee2(2,:) =  r_ee2(2,:)/dble(n_chain/2)
              r_ee2(3,:) =  r_ee2(3,:)/dble(n_chain_d)
              
              r_ee_mean(:,:) = r_ee_mean(:,:) + r_ee(:,:)
              r_ee_mean2(:,:) = r_ee_mean2(:,:) + r_ee2(:,:)
!

        case(3)    
              r_ee_mean(:,:) = r_ee_mean(:,:) / real(n_count)
              r_ee_mean2(:,:) = r_ee_mean2(:,:) / real(n_count)

              print '(/15x,a/)', "End to End Radius  "
              print '(a30,3f16.7)',"Top brushes, Ree= ",r_ee_mean(1,:)
              print '(a30,3f16.7)',"Bottom brushes, Ree= ",r_ee_mean(2,:)
              print '(a30,3f16.7)',"Melt,  Ree= ",r_ee_mean(3,:)
! ------     Square end to end distance:
              print '(/a30,4f16.7)',"Top brushes, Ree2= ",r_ee_mean2(1,:)   ,sum(r_ee_mean2(1,:))
              print '(a30,4f16.7)',"Bottom brushes, Ree2= ",r_ee_mean2(2,:) ,sum(r_ee_mean2(2,:)) 
              print '(a30,4f16.7/)',"Melt,  Ree2= ",r_ee_mean2(3,:)         ,sum(r_ee_mean2(3,:))
         end select             

        end subroutine r_end_to_end_2v
!
!                                **** End to End Radius  ****
!
        subroutine r_end_to_end(mode,i_case,n_mon,spabs,r_ee,r_ee_mean,r_ee2,r_ee_mean2)   
        integer , intent(in) :: mode,n_mon,i_case
        real (kind=8) , intent(in)   :: spabs(:,:)
        real (kind=8) , intent(out)  :: r_ee(:,:),r_ee_mean(:,:),r_ee2(:,:),r_ee_mean2(:,:)
        integer,  save :: n_count(3),n_case
        real (kind=8)  :: v_dummy(3) , r_dummy(3),dist
        integer :: i_chain, n_chain,i_part,i

        select case(mode)
        case(1)  !  init
            n_case = 3 ! size(r_ee,dim=1)
!          allocate( r_ee_mean(n_case,3) ,n_count(n_case) ) 
          n_count(:) = 0
          r_ee_mean(:,:) = 0.
          r_ee_mean2(:,:) = 0.
          open(unit=46,file='R_ee_dist.mide',status='unknown')
          write(46,*)  '#Distance of Ree for brushes vs time'
        case(2)  !  calculation
              n_count(i_case) = n_count(i_case) + 1
              n_chain = size (spabs,dim = 1) /n_mon 
!              print*,"n_chain=",n_chain
              r_ee(i_case,:) = 0.d0
              r_ee2(i_case,:) = 0.d0
!
!        **** R end to end = Ree **** 
!
              do i_chain = 1,n_chain
                 v_dummy(:) = spabs(i_chain*n_mon, :) - spabs((i_chain-1)*n_mon+1,:) 
                  r_ee(i_case,:)=r_ee(i_case,:)+ v_dummy(:) 
                  r_ee2(i_case,:)=r_ee2(i_case,:)+ v_dummy(:)**2 
!!!! WARN: commented for speed-up
!                  r_ee_vec(i_chain,:) = v_dummy(:)
!!!             if(i_case==1.or.i_case==2) then  ! for brushes
!!!                 dist = sqrt( dot_product( v_dummy(:),v_dummy(:) ) ) 
!!!                 write(46,*) n_count(i_case),dist
!!!             end if
              end do
              r_ee(i_case,:) =  r_ee(i_case,:)/dble(n_chain)
              r_ee2(i_case,:) =  r_ee2(i_case,:)/dble(n_chain)
              r_ee_mean(i_case,:) = r_ee_mean(i_case,:) + r_ee(i_case,:)
              r_ee_mean2(i_case,:) = r_ee_mean2(i_case,:) + r_ee2(i_case,:)
!                
!

        case(3)    
              do i = 1,n_case
              r_ee_mean(i,:) = r_ee_mean(i,:) / real(n_count(i))
              r_ee_mean2(i,:) = r_ee_mean2(i,:) / real(n_count(i))
              end do
              print '(/15x,a/)', "End to End Radius:  "
              print '(a30,3f16.7)',"Top brushes, Ree= ",r_ee_mean(1,:)
              print '(a30,3f16.7)',"Bottom brushes, Ree= ",r_ee_mean(2,:)
              print '(a30,3f16.7)',"Melt,  Ree= ",r_ee_mean(3,:)
! ------     Square end to end distance:
              print '(/a30,4f16.7)',"Top brushes, Ree2= ",r_ee_mean2(1,:)  ,sum(r_ee_mean2(1,:))
              print '(a30,4f16.7)',"Bottom brushes, Ree2= ",r_ee_mean2(2,:),sum(r_ee_mean2(2,:))
              print '(a30,4f16.7/)',"Melt,  Ree2= ",r_ee_mean2(3,:)        ,sum(r_ee_mean2(3,:))
         end select             

        end subroutine r_end_to_end
!        
!  *****         Angle with a plane        
!
        function  angle_plane(vector,p_normal) result(angle)
        real (kind=8) :: vector(3), p_normal(3),angle 
        real (kind=8)  :: pi 
        pi = 4.*atan(1.)
         
        angle = pi/2. - acos ( dot_product(p_normal,vector)/         &
        sqrt ( dot_product(vector,vector)*dot_product(p_normal,p_normal) ) )
        angle = angle* 180./pi
        end function angle_plane
!
!   ***** Drop to a file the velocity history of the last brush monomer
!
  subroutine drop_v_end(mode,r_time,n_chain,n_mon)
  integer ,intent(in) :: mode,n_chain,n_mon 
  real(kind=8) ,intent(in) :: r_time
  integer :: i,j
  select case (mode)
  case(1)
      open(unit=47,file='v_end_history.mide',status='unknown')
      write(47,*) '# Velocity of the last monomer in the brush ' 
  case(2)
      do i=1,n_chain
      j = i*n_mon 
      write(47,'(4f16.7)') r_time,v(j,:) ! Be careful because top and botom walls vel. have diff. signs
      end do
  end select
  end subroutine drop_v_end
!
! ** Bond orientation calculations
!
  subroutine bond_orientation(mode,D,n_mon,n_chain,n_mon_d,n_chain_d,r0_unfold,count_bonds,bond_orient_histo)
  real (kind=8), intent(in) :: D
  integer      , intent(in) :: n_mon,n_chain,n_mon_d,n_chain_d,mode 
  real (kind=8), intent(in) :: r0_unfold(:,:)
  real (kind=8),intent(inout) :: bond_orient_histo(:,:)
  integer, save :: n_layers,n_chains,l_count,melt0
  integer,  parameter :: pi= 3.141516
  integer       :: i_chain,i_mon,i_layer,i,at1,j
!  integer, allocatable:: count_bonds(:,:)
!  integer   :: count_bonds(2,size(bond_orient_histo(:,:),dim=2))= 0 
  integer   :: count_bonds(:,:) 
  real (kind=8), save :: dz
  real (kind=8)       :: bond(3),z1,length_2,s1,s2
  select case (mode)
  
  case(1)
      bond_orient_histo(:,:) = 0.0
      n_layers = size(bond_orient_histo(:,:),dim=2)
!debug      print *,n_layers ; stop
!      allocate( count_bonds(2,n_layers))
!      n_chains = size(r0,dim=1)/n_mon
      dz = D/real(n_layers)
      open(unit=49,file='bond_orient.mide',status='unknown')
      write(49,*)  '# z      0.5*(3*<cos(theta)^2>-1) angle of bond with the 001 vector'
      l_count = 0
      melt0 = n_mon*n_chain
      count_bonds(:,:) = 0 
  case(2)
!not used by now      l_count = l_count + 1  ! n times we are here
      
!      print*,count_bonds
!      stop
! Brush      
       do i_chain =1,n_chain
       do i_mon= 1, n_mon-1

          at1 = (i_chain-1)*n_mon+i_mon
!          print *,"at1=",at1,at1+1
          z1= (r0_unfold(3,at1+1) + r0_unfold(3,at1))/2.
       do i_layer = 1 , n_layers
       if (z1<real(i_layer)*dz .and. z1>=real(i_layer-1)*dz) then
           count_bonds(1,i_layer) = count_bonds(1,i_layer) + 1
           bond(:) = r0_unfold(:,at1+1) - r0_unfold(at1 ,:)
           length_2 = dot_product(bond,bond)
!           print*,'length_2=',length_2
!           if(length_2>9.) then
!               print *,"bond too high"
!               stop
!            endif
!           bond(:) = bond(:) / sqrt( dot_product(bond,bond) )
           bond_orient_histo(1,i_layer) = bond_orient_histo(1,i_layer) +     & 
                                   bond(3)**2/length_2   ! melt
!           print*,"cos^2(theta)=" ,(bond(3))**2/length_2                        
           exit  
       end if

       end do
       end do
       end do
! Melt 
       do i_chain =1,n_chain_d
       do i_mon= 1, n_mon_d-1

           at1 = melt0 + (i_chain-1)*n_mon+i_mon
!           print*,at1
           z1= (r0_unfold(3,at1+1) + r0_unfold(3,at1))/2.
       do i_layer = 1 , n_layers
       
       if (z1<real(i_layer)*dz .and. z1>=real(i_layer-1)*dz) then
           count_bonds(2,i_layer) = count_bonds(2,i_layer) + 1
           bond(1:3) = r0_unfold(at1 + 1,:) - r0_unfold(:,at1)
           length_2 =  dot_product(bond,bond) 
           bond_orient_histo(2,i_layer) = bond_orient_histo(2,i_layer) + & 
                                 bond(3)**2/length_2   ! melt

           exit  
       end if

       end do
       end do
       end do
       
  case(3)
      ! this is gramatically correct
!       where(count_bonds /= 0) bond_orient_histo(:,:) =  bond_orient_histo(:,:) / real(count_bonds(:,:)) 
!     bond_orient_histo(:,:) =  bond_orient_histo(:,:) /real(l_count ) ! redundant with the other normaliz ? YES
print '(/a/)'," * bond orientation note: if not counts in a given layer, S is set to 0  " 
     do i = 1,n_layers
     if(bond_orient_histo(1,i) == 0.) then
         s1 = 1./3.
     else 
         s1 = bond_orient_histo(1,i) / real(count_bonds(1,i))
     end if
     if(bond_orient_histo(2,i) == 0.) then
         s2 = 1./3.
     else 
         s2 = bond_orient_histo(2,i)/ real(count_bonds(2,i))
     end if
! write file     
     write(49,'(3f16.7)') dz/2. + real(i-1)*dz,       &
     0.5*(3.*s1-1.),0.5*(3.*s2 - 1.)
! set up final bond_orient_histo:     
     bond_orient_histo(1,i) = 0.5*(3.*s1-1.) 
     bond_orient_histo(2,i) = 0.5*(3.*s2-1.) 



     end do
  end select     
  end subroutine bond_orientation
!
! Profile calculation routine of order parameter 
!
   subroutine profile_orient(D,vec_in,r_cms,prof,mean_prof,l_count)
   
   real (kind=8), intent(in) :: vec_in(:,:),r_cms(:,:)
   real (kind=8), intent(in) :: D
   integer, intent(inout) :: l_count(:) 
   real (kind=8), intent(out) :: prof(:),mean_prof(:)
   real (kind=8) :: dz,length_2
   integer :: i_vec,j,k,n_layers,i_layer,vec_dim

   n_layers = size(prof)
   dz = D/real(n_layers)
   vec_dim = size(vec_in,dim=1)
   prof(:) = 0.
   do i_vec = 1, vec_dim
    do i_layer = 1, n_layers
    if (r_cms(i_vec,3)<real(i_layer)*dz .and. r_cms(i_vec,3) > real(i_layer-1)*dz ) then
        l_count(i_layer) = l_count(i_layer) + 1
        
 ! calculation of the order parameter        

           length_2 = dot_product(vec_in(i_vec,:),vec_in(i_vec,:) ) 
! I use now this: Rz^2-(Rx^2+Ry^2)/2 / length^2
        prof(i_layer) = prof(i_layer) +         &
        ( vec_in(i_vec,3)**2 - 0.5*(vec_in(i_vec,1)**2 + vec_in(i_vec,2)**2) ) / length_2   
        exit
    end if
    end do
   end do
       mean_prof(:) = mean_prof(:) + prof(:)
   
   end subroutine profile_orient
   
!
! Writes a file of the vector and the name of the argument


   subroutine write_file(f_name,comment,x_axes,y_axes)
   character (len=*), intent(in) :: f_name, comment
   real(kind =8) ,intent(in) :: x_axes(:),y_axes(:)
   integer :: i 
   open(unit=92,file=f_name,status='unknown')
   write(92,*) comment
   do i=1,size(y_axes,dim=1) 
   write(92,'(2f16.7)')  x_axes(i),y_axes(i)
   end do
   close(92) 
   end subroutine write_file

!
! CM histograms 
!
subroutine histo_cm(mode,D,n_chain,n_chain_d,vector,histo_out)

implicit none
real (kind=8), dimension(:,:), intent(in) :: vector
real (kind=8) , intent(in) :: D 
integer , intent(in) :: mode,n_chain,n_chain_d 
real (kind=8) , dimension(:,:), intent(out) :: histo_out 
real (kind=8),save  :: dz
integer, save :: l_count, n_part, n_dim
integer :: n_max,i_step,i_part,n_step_clau

select case (mode)
case(1)
    histo_out(:,:) = 0.
    n_dim = size(histo_out,dim=2)
    dz = D/real(n_dim)
    n_part = size(vector,dim=1) 
    l_count = 0
case(2)
    l_count = l_count + 1
! Brushes
    do i_part = 1, n_chain 
        do i_step = 1, n_dim
        
        if(vector(i_part,3) < dz*i_step .and. vector(i_part,3) >= dz*(i_step-1) ) then
          histo_out(1,i_step) = histo_out(1,i_step) + 1
          exit
        end if
        end do
    end do
!  Melt
    do i_part = n_chain, n_chain+n_chain_d 
        do i_step = 1, n_dim
        
        if(vector(i_part,3) < dz*i_step .and. vector(i_part,3) >= dz*(i_step-1) ) then
          histo_out(2,i_step) = histo_out(2,i_step) + 1
          exit
        end if
        end do
    end do
    ! note: the normalization and output is done in the main program     
case(3)
     histo_out(:,:) = histo_out(:,:) /real (l_count ) ! here is normalized to the number of time steps
    
end select
end subroutine histo_cm

! Calculates the mean value of Ree or R in the middle of the sample

subroutine mean_R_center_box(z,z_int,R_vec,R_cm,mean_R,count_R)

real (kind=8) , dimension(:,:), intent(in) :: R_vec,R_cm
real (kind=8) , intent(in) :: z,z_int
integer, intent(inout) :: count_R
real(kind=8) ,intent(inout) :: mean_R(3)
integer :: i,n_part

!debug print*,r_vec(:,:) ; stop
    
n_part = size(R_vec(:,:),dim=1)

    do i = 1 , n_part
    if (R_cm(i,3) < z/2.+z_int/.2 .and. R_cm(i,3) > z/2.-z_int/.2 ) then
        
        mean_R(:) = mean_R(:) + R_vec(i,:)
        count_R = count_R+1

    end if
    end do
    

end subroutine mean_R_center_box


subroutine calc_cms(mode) !,spabs_l,r_cm_vec)
        integer , intent(in) :: mode
!        real (kind=8) ,intent(in) :: spabs_l(:,:)
!        real (kind=8) ,intent(out) :: r_cm_vec(:,:)
        integer,  save :: n_count
        real (kind=8)  :: v_dummy(3) , r_cm_local(3), r_dummy(3)
        integer :: i_chain, i_part, i, j
        integer, save :: n_case,n_free_chains
        real (kind=8), save :: inv_n_mon,inv_n_mon_d

      select case(mode)

        case(1)  !  init

            n_free_chains = n_chain_d ! warning: I do not use anymore explicity vapor chains + n_chain_vap
            inv_n_mon = 1./dble(n_mon)
            inv_n_mon_d = 1./dble(n_mon_d)

        return

        case(2)  !  calculation

!  CM calculation

             do i_chain = 1,n_chain                    ! brush
!old    
              r_cm_local(:) = sum(r0_unfold(:, (i_chain-1)*n_mon+1:i_chain*n_mon),dim=2)*inv_n_mon  
                 
                  r_cms(i_chain,:) = r_cm_local(:)*inv_n_mon
!                print*,r_cms(i_chain,:)
             end do
             do i_chain = n_chain+1, n_chain+n_free_chains ! melt/droplet
!old
                r_cm_local(:) =  &
                sum(r0_unfold(:, (i_chain-1)*n_mon_d+1:i_chain*n_mon_d),dim=2) *inv_n_mon_d

                  r_cms(i_chain,:) = r_cm_local(:)*inv_n_mon_d
             end do
         end select

        end subroutine calc_cms



#   endif !  PROFILES 

#   if SYSTEM==1 
!
!  ROUTINES FOR DROPLET SIMS
!
!
!  * Rebuilds the droplet using info from rebuild2 routine and the Rg_x
!  * Also writes out Rg_x and _z versus time 
!
! mode = 0=init. /=0 box position of the CM
subroutine rebuild_drop(mode,rg_x) !,r_coors,r_cm )
integer, intent(in)  :: mode 
real (kind=8) , intent(out) :: rg_x
real (kind=8) :: bound 
real (kind=8), save  :: shift_cm(3),box_center
real (kind=8), allocatable ,save  :: r0_aux(:,:)
integer :: i_cm,int_curr,int_next,i
real (kind=8) :: rg_x_ori
integer, save :: icount, n_cha
logical, parameter :: debug=.false.

select case (mode) 

case(0)  ! INIT
    
! For the first step: the guessed CM is used. Now is calculated in
! rebuild_drop2. 
! NOTE: only the position in X is used 
 
!  shift_cm(:) =(/ guessed_drop_cm(1), 0._8 , 0._8  /)
!  print '(a,3f10.5/)','[rebuild_drop]  *  First drop CM =  ',shift_cm
        icount = 0 
! Number of chains         
        n_cha = part_init_d/n_mon + nm/n_mon_d
        box_center = boundary(1)/2.
        allocate( r0_aux(n_mon_tot,3) )

    case(1)   !   CM in the left third of the box 
        icount = icount + 1
! * The algorithm assumes that the droplet smaller than half the box in X. 
! * NOTE: this translates brush and free chains     

!  ---- Calculates Rg_x before translation

           rg_x_ori = 0. 
           do i = part_init_d+1,n_mon_tot !melt
               rg_x_ori = rg_x_ori + ( r0_unfold(1,i) - drop_cm(1) )**2
           end do
           rg_x_ori = rg_x_ori/dble(nm)
     if (debug) print *,"Rg ori = ", rg_x_ori

!  ---- Translate           
!Modified by Kevo
           do i_cm = 1, n_chain ! loop over brush
                 int_curr = (i_cm-1)*n_mon  + 1 
                 int_next = (i_cm)*n_mon   
             if (  r_cms(i_cm,1) > box_center    ) then
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) - boundary(1)
             else
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) 
             end if
           end do

          do i_cm = 1, n_chain_d ! loop over melt
                 int_curr = part_init_d + (i_cm-1)*n_mon_d  + 1 
                 int_next = part_init_d + (i_cm)*n_mon_d   
             if (  r_cms(i_cm,1) > box_center    ) then
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) - boundary(1)
             else
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) 
             end if
           end do
!  ---- Calculates Rg_x after translation

           rg_x = 0. 
           do i = part_init_d+1,n_mon_tot
               rg_x = rg_x + ( r0_aux(i,1) - drop_cm(1) )**2
           end do
           rg_x = rg_x/dble(nm)
     if (debug) print *,"Rg new  = ", rg_x
!  ---- Decide if translation is done or not 

           if( rg_x < rg_x_ori ) then ! translation is 
              r0_unfold(1,:) = r0_aux(:,1)
           end if 
           ! else translation is not done
     case(2)   !   CM in the central third of the box 
         ! nothing to do
           rg_x = 0. 
           do i = part_init_d+1,n_mon_tot
               rg_x = rg_x + ( r0_unfold(1,i) - drop_cm(1) )**2
           end do
           rg_x = rg_x/dble(n_mon_tot-part_init_d) ! MORE EFFICIENT LATER
           return
     case(3)   !   CM in the right third of the box 

!  ---- Calculates Rg_x before translation

           rg_x_ori = 0. 
           do i = part_init_d+1,n_mon_tot
               rg_x_ori = rg_x_ori + ( r0_unfold(1,i) - drop_cm(1) )**2
           end do
           rg_x_ori = rg_x_ori/dble(nm)
     if (debug) print *,"Rg ori = ", rg_x_ori
!  ---- Translate           

           do i_cm = 1, n_chain !loop over brush
                 int_curr = (i_cm-1)*n_mon  + 1 ! monomero inicial dentro de la cadena i_cm
                 int_next = (i_cm)*n_mon   ! monomero final dentro de la cadena i_cm

             if (  r_cms(i_cm,1) < box_center    ) then
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) + boundary(1)
             else
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) 
             end if
           end do

           do i_cm = 1, n_chain_d !loop over melt
                 int_curr = part_init_d + (i_cm-1) * n_mon_d  + 1 ! monomero inicial dentro de la cadena i_cm
                 int_next = part_init_d + (i_cm) * n_mon_d     ! monomero final dentro de la cadena i_cm
                 if (debug) print *,"int_curr = ", int_curr, "int_next = ", int_next

             if (  r_cms(i_cm,1) < box_center    ) then
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) + boundary(1)
             else
                 r0_aux(int_curr:int_next,1) =  r0_unfold(1,int_curr:int_next) 
             end if
           end do


!  ---- Calculates Rg_x after translation

           rg_x = 0. 
           do i = part_init_d+1,n_mon_tot
               rg_x = rg_x + ( r0_aux(i,1) - drop_cm(1) )**2
           end do
           rg_x = rg_x/dble(nm)
     if (debug) print *,"Rg new = ", rg_x
!  ---- Decide if translation is done or not 

           if( rg_x < rg_x_ori ) then ! translation is done 
              r0_unfold(1,:) = r0_aux(:,1)
              if(debug) print*,"Translation is done"
           else ! as it is not done, we store Rg_x as the ori value 
               rg_x = rg_x_ori 
           end if 
           ! else translation is not done

! Substract Lx/2 in order to have the droplet CM in the first half of the BOX.
! (needed for the future translations)

!                r0_unfold(1,:) = r0_unfold(1,:) - box_center

      case default
             print *,'[rebuild drop] Incorrect case. mode value incorrect'
             stop
      end select


end subroutine rebuild_drop

subroutine calc_drop_rg(rg_x)
        real(kind=8),intent(in) :: rg_x
        real(kind=8) :: rg_z
        integer :: i

! Calculates Rg_z and writes out Rg vs. time 

           rg_z = 0. 
           do i = part_init_d+1,n_mon_tot
           rg_z = rg_z + ( r0_unfold(3,i) - drop_cm(3) )**2
           end do
           rg_z = rg_z/dble(n_mon_tot-part_init_d)

! Writes out            
           write(55,'(i7,x,2(g17.10,x))') i_time,rg_x,rg_z 

end subroutine calc_drop_rg


subroutine rebuild_drop2(mode)

integer, intent(in) ::  mode
real (kind=8) , allocatable, save :: scaled_3_x(:),scaled_4_x(:),r_local(:,:)
integer :: n3_box(3), n4_box(4),i, index_3,index_4
integer,save :: cm_ok
logical, parameter :: debug=.false.
!below, variables for case(2). More accuarate drop_cm calculation.
logical :: control
real (kind=8) :: x_cm, x_cm_new, delta_cm, delta_cm_new
real (kind=8),dimension(nm) :: r0_cent
!below, variables for case(3). Optimized drop_cm calculation
integer :: n_boxes, i_index, max_box
real (kind=8) :: delta_L = 0.25
integer, allocatable, save :: box_count(:)
!below, variables for case(4). Bisection drop_cm calculation 
integer :: case1, case2
real (kind=8) :: box_min, box_middle, box_max

! Find in which third of the box there are more beads
! N_CM take the values 1,2,or 3 according with the third in the box that has
! more beads  
select case (mode)

case(0)  ! init 
    n_boxes = int( boundary(1) / delta_L ) + 1
    allocate(box_count(n_boxes))
    allocate( scaled_3_x(nm),  scaled_4_x(nm) , r_local(n_dim,n_mon_tot))
!Corrected by Kevo 4/2016

    n3_box(:) = 0 
    n4_box(:) = 0 
    r_local(:,:) = 0

     print '(/a/)' ,'  * Finding first droplet CM  ' 

case(1) ! Calculation for the first time 

! In which third or fourth of the box there are mor beads ?     

    n3_box(:) = 0 
    n4_box(:) = 0 

    scaled_3_x(:) = 3._8*r0(1,part_init_d+1:n_mon_tot) * inv_boundary(1)
    scaled_4_x(:) = 4._8*r0(1,part_init_d+1:n_mon_tot) * inv_boundary(1)
        do i = 1,nm
         index_3 = int( scaled_3_x(i) )  + 1
         index_4 = int(scaled_4_x(i))  + 1


         n3_box(index_3) =  n3_box(index_3) + 1 
         n4_box(index_4) =  n4_box(index_4) + 1 

        end do 
       cm_ok = maxloc(n3_box(:),dim=1)
       
! Global, used in  rebuild_drop  

       n_cm = cm_ok

   if(debug) then
       print '(a,3i)' ,"        Number of bead per box (3 boxes): ",  n3_box
       print '(a,4i)' ,"        Number of bead per box (4 boxes): ",  n4_box 
       print '(a,i/)' ,"        DROP CM should be in box ", cm_ok 
   end if

! The brush goes always in the same position 

        r_local(:,:)  = r0(:,:)
! Rebuild the drop 
     select case (cm_ok) 

     case(1) ! Most of the droplet in the first box | o |   |   | 

        do i = 1,nm
         index_3 = int( scaled_3_x(i) )  + 1
         index_4 = int(scaled_4_x(i))  + 1
             if(index_4 == 3 .or. index_4 == 4) then
                 r_local(1, part_init_d+i ) = r_local(1, part_init_d+i) - boundary(1)
             end if
        end do
     case(2) ! Most of the droplet in the second box |   | o |   | 

            r_local(1, part_init_d+1:n_mon_tot) =  r0(1,part_init_d+1:n_mon_tot)
     case(3) ! Most of the droplet in the third box  |   |   | o | 

            do i = 1,nm
                index_3 = int( scaled_3_x(i) )  + 1
                index_4 = int(scaled_4_x(i))  + 1
                if(index_4 == 1 .or. index_4 == 2) then
                    r_local(1, part_init_d+i) = r_local(1, part_init_d+i) + boundary(1)
                end if 
            end do
     end select

! Calculate the CENTER OF MASS   (global var.)   

        drop_cm(:) = sum(r_local(:,part_init_d+1:n_mon_tot),dim=2)/dble(nm)

!        print '(a,4f10.5/)',"[rebuild_drop2] Estimated droplet CM: ",drop_cm(:)

case(2) !More accuarate drop_cm calculation, for systems 
        ! with a non negligible vapor fraction. Works for 1 droplet only. 
        !Also gives r0_cent, wich is r0 with the center of mass inx in
        !boundary(1)/2. by Kevo 5/2016.
        !The idea is to calculate the center of mass of the melt. Then translate
        !the system, so that the centre of mass is in Lx/2. Then apply PBC to
        !have all particles in 0 < x < Lx. Then calculate the center of mass of
        !the melt. Then translate the system, so that the centre of mass is in
        !Lx/2 and so on... This process has a fixed point, when the centre of
        !mass of the melt is in Lx/2, which for 1 droplet, coincides with the
        !center of mass of the droplet. WARNING: arbitrary convergence criteria

    !Define r0_cent for the first time in this step:
    r0_cent(:) = r0(1,part_init_d+1:n_mon_tot)
    x_cm = sum(r0_cent(:))/dble(nm) 
    control = .true.
    drop_cm(1) = x_cm
!While not converged
    do while(control) 
        r0_cent(:) = r0_cent(:) - x_cm + boundary(1)/2. !Translate whole system so that 
                                                            ! x_cm is in Lx / 2  

        !Apply PBC for droplet only
        do i_part = 1,nm
            if (r0_cent(i_part).lt.0.) then
                r0_cent(i_part) = r0_cent( i_part) + boundary(1)
            else if (r0_cent(i_part).gt.boundary(1)) then
                r0_cent( i_part) = r0_cent( i_part) - boundary(1)
            end if
        end do    
    
        !get new centre of mass x coordinate after aplying PBC
        x_cm_new = sum(r0_cent(:))/dble(nm)
         
        delta_cm =      x_cm -  boundary(1) / 2.      !Quantify translation 
        delta_cm_new =  x_cm_new - boundary(1) / 2.  
        
        !Check if converged
        if(abs(delta_cm_new).le.abs(delta_cm)) then !Transaltion in x mus be less than previous transaltion
            if (abs(delta_cm_new).le.0.25) then !The distance between centre of mass
                control = .false.        !and box centre shall be less than the numer in this line
            end if                      ! HARDCODE WARNING: 0.25 in sigma units            
        end if     
        x_cm = x_cm_new 
        drop_cm(1) = drop_cm(1) +  delta_cm_new                
    end do !End while
    
    drop_cm(2:3) = sum(r0(2:3,part_init_d+1:n_mon_tot),dim=2)/dble(nm)
    n_cm = int(drop_cm(1) * inv_boundary(1) * 3.) + 1 

case(3) !Another way of calculating drop_cm. The idea is to divide Lx in many
        !boxes, and count the number of melt particles in each box. Then
        !transalte the system so that the center of the droplet is at Lx/2
    delta_L = 0.25
    n_boxes = boundary(1)/delta_L
    box_count(:) = 0
    do i_part = part_init_d+1, n_mon_tot !loop over all melt particles 
        i_index = int( r0(1,i_part) * inv_boundary(1) * n_boxes + 1 )
        box_count(i_index) = box_count(i_index) + 1
    end do
    max_box = maxloc(box_count,dim=1)
    drop_cm(1) = ( dble(max_box) - 0.5 ) * delta_L 

    drop_cm(2:3) = sum(r0(2:3,part_init_d+1:n_mon_tot),dim=2)/dble(nm)
    n_cm = int(drop_cm(1) * inv_boundary(1) * 3.) + 1

case(4) !Yet another way to calculate drop_cm. This case is by bisection.
    box_min = 0. !Define first box min and max
    box_max = boundary(1)
    do while((box_max-box_min).gt.0.25) !WARNING: Resolution to drop_cm is HARDCODED
        box_middle = (box_min + box_max) / 2.
        case1 = 0
        case2 = 0
        do i_part = part_init_d+1, n_mon_tot !loop over all melt particles
            if((r0(1,i_part).ge.box_min).and.(r0(1,i_part).lt.box_middle)) then! Particle is on the left side of the box
                case1 = case1 + 1
            else if ((r0(1,i_part).ge.box_middle).and.(r0(1,i_part).lt.box_max)) then
                                 ! Particle is on the right side of the box
                case2 = case2 + 1
            end if
        end do        
        if(case1.gt.case2) then !More particles on the middle left side of the box
            box_max = box_middle
        else if (case2.gt.case1) then! More particles on the middle right side of the box
            box_min = box_middle
        else if(case1.eq.case2) then !Equal numer of particles on the right and on the left
            box_max = box_middle
            box_min = box_middle
        end if
    end do    
    drop_cm(1) = box_middle
    drop_cm(2:3) = sum(r0(2:3,part_init_d+1:n_mon_tot),dim=2)/dble(nm)
    n_cm = int(drop_cm(1) * inv_boundary(1) * 3.) + 1

!OUTPUT : drop_cm and n_cm
  end select ! mode

end subroutine rebuild_drop2
!
! Dens. Profile:  for the droplet and brush layer: 2D profile
! 
subroutine dens_prof(mode,histo,histo_b,histo2,histo_b2,n_box,r_box)

     integer, intent(in) :: mode,n_box(3)
     real (kind=8), dimension(:,:,:) ,intent(inout) :: histo, histo_b, histo2,histo_b2
     real (kind=8), intent(in) :: r_box(3)
     real (kind=8) :: r_box_max(3),r_box_min(3)
     real (kind=8), allocatable, save :: r_scaled(:,:), histo_step(:,:),histo_step_b(:,:)
     real (kind=8), save :: inv_r_box(3)
!     integer, allocatable,save :: histo_f(:)
!     real (kind=8)       :: histo_step(n_box(1),n_box(3)),histo_step_b(n_box(1),n_box(3))
     integer :: i,j,k,i_part,r_index(3)
     integer, save :: i_time,nx,ny,nz
     integer :: f_index_x
     integer, parameter :: debug=.false.

select case (mode)
case(0) ! Init.
    print "(/a/)", "  *   Initialising dens_prof"
    histo(:,:,:) = 0.0
    histo_b(:,:,:) = 0.0

    histo2(:,:,:) = 0.0
    histo_b2(:,:,:) = 0.0

    nx = size(histo,dim=1)
    ny = size(histo,dim=2)
    nz = size(histo,dim=3)

    inv_r_box(:) = 1./r_box(:)

    i_time = 0 
    print *, "  *  nx,ny,nz,n_part   = ",nx,ny,nz,n_part

#ifdef FORCE_PROFILE
    force_prof(:,:) = 0.0
    allocate ( histo_f(nx))
    histo_f(:) = 0
#endif

    allocate (r_scaled(3,n_part))
    allocate (  histo_step(nx,nz),histo_step_b(nx,nz) )

    r_scaled(1:3,:) = 0.0

case(1) ! normalization and writing


!    print*,histo(:,:,:)

    histo(:,:,:) = histo (:,:,:) / (dble(i_time)* r_box(1)*r_box(2)*r_box(3) ) 
    histo_b(:,:,:) = histo_b (:,:,:) / (dble(i_time)*  r_box(1)*r_box(2)*r_box(3) ) 
    
    histo2(:,:,:) = histo2 (:,:,:) / (dble(i_time)*  r_box(1)*r_box(2)*r_box(3) ) 
    histo_b2(:,:,:) = histo_b2 (:,:,:) / (dble(i_time)*  r_box(1)*r_box(2)*r_box(3) ) 


      print '(/a)', "  *  DROP: Writing num_density.mide ... " 
      print '(a/)', "  *  BRUSH:Writing num_density_b.mide ... " 

!---  Write out droplet profile 

    open(unit=35,file="num_density.mide",status="unknown")
! Header for automatic reading from tecplot !OLD not in use anymore
        write(35,*) 'variables= "x","z","density"'
        write(35,*) 'zone i= ',nx,' ,j= ',nz,' ,f=point'
        write(35,*) 

! droplet
        do k = 1,nz
            do i = 1,nx
                write(35,'(4f16.8)')    dble(i-1)*r_box(1) + r_box(1)/2.,                        &
                    dble(k-1)*r_box(3) + r_box(3)/2., histo(i,1,k) 
            enddo
        enddo
    close(35)

!---  Write out dispersion of droplet profile  sqrt( <dens2>-<dens>^2 )

    open(unit=35,file="num_sig_density.mide",status="unknown")
! Header for automatic reading from tecplot
        write(35,*) 'variables= "x","z","sig(dens)"'
        write(35,*) 'zone i= ',nx,' ,j= ',nz,' ,f=point'
        write(35,*) 

! droplet sigma
        do k = 1,nz
                  do i = 1,nx
!         do j = 1,ny
                write(35,'(4f16.8)')    dble(i-1)*r_box(1) + r_box(1)/2.,                        &
                                        dble(k-1)*r_box(3) + r_box(3)/2.,                        &
                                       sqrt( histo2(i,1,k)- histo(i,1,k)**2 ) 
                enddo
!        enddo
    enddo
    close(35)

!---  Write out brush profile 

    open(unit=35,file="num_density_b.mide",status="unknown")
! Header for automatic reading from tecplot
        write(35,*) 'variables= "x","z","density"'
        write(35,*) 'zone i= ',nx,' ,j= ',nz,' ,f=point'
        write(35,*) 

! droplet
        do k = 1,nz
             do i = 1,nx
!         do j = 1,ny

                write(35,'(4f16.8)')    dble(i-1)*r_box(1) + r_box(1)/2.,  &
                                        dble(k-1)*r_box(3) + r_box(3)/2., histo_b(i,1,k) 
                enddo
!        enddo
    enddo
    close(35)

!---  Write out brush  dispersion profile:  sqrt( <dens2>-<dens>^2 )

    open(unit=35,file="num_sig_density_b.mide",status="unknown")
! Header for automatic reading from tecplot
        write(35,*) 'variables= "x","z","density"'
        write(35,*) 'zone i= ',nx,' ,j= ',nz,' ,f=point'
        write(35,*) 

! Droplet

      do k = 1,nz
!         do j = 1,ny
          do i = 1,nx
                write(35,'(4f16.8)')    dble(i-1)*r_box(1) + r_box(1)/2.,  &
                                        dble(k-1)*r_box(3) + r_box(3)/2.,  &
                                        sqrt( histo_b2(i,1,k) - histo_b(i,1,k)**2 )
          enddo
!        enddo
      enddo
    close(35)

#ifdef FORCE_PROFILE    
!
! Write out force profile in X 
!
        do i = 1,nx
            if(histo_f(i)> 0) then
                force_prof(i,:) = force_prof(i,:)/dble(histo_f(i)*i_time)
            else
                force_prof(i,:) = 0.
            end if
        end do
    
    open(unit=35,file="force_profile.mide",status="unknown")
    write(35,*) "#x_box   Fx       Fy     Fz   [on brush heads ]"
    do i = 1,nx
        write(35,'(4f16.8)')    dble(i-1)*r_box(1),force_prof(i,:) 
    end do
    close(35)

#endif

case(2) ! calculates profile droplet and brush ! [CURRENTLY NOT USED see case 3]
! time counter
      i_time = i_time + 1

      ! VERY UNEFFICIENT(but it works):   check order N algorithms case=3 
! 
        do i = 1,nx
            r_box_min(1) = r_box(1)*real(i-1)
            r_box_max(1) = r_box(1)*real(i)
            do k = 1,nz
                r_box_min(3) = r_box(3)*real(k-1)
                r_box_max(3) = r_box(3)*real(k)
                do i_part = 1,n_part
                    if ( r0_unfold(1,i_part)>=r_box_min(1) .and. r0_unfold(1,i_part)<r_box_max(1) .and.   &
                        r0_unfold(3,i_part)>=r_box_min(3) .and. r0_unfold(3,i_part)<r_box_max(3)  ) then
                    ! droplet                             
                    if(a_type(i_part)==3) then  ! droplet 
                        histo(i,1,k) =  histo(i,1,k) + 1 
                    end if
                    ! brush                       
                    if(a_type(i_part)==2.or.a_type(i_part)==1) then ! brush grafted head and tail
                        histo_b(i,1,k) =  histo_b(i,1,k) + 1 
                    end if
                end if
            enddo  ! particles
        enddo
    enddo

case(3) ! trial of order N algorithm  ! WORKS ! 
! time counter
      i_time = i_time + 1
      histo_step(:,:) = 0.
      histo_step_b(:,:) = 0.

      ! Change scale so that int() means find the right binning box for the particle      

      r_scaled(1,:) = r0_unfold(1,:)*inv_r_box(1)
      r_scaled(2,:) = r0_unfold(2,:)*inv_r_box(2)
      r_scaled(3,:) = r0_unfold(3,:)*inv_r_box(3)

      do i_part = 1,n_part

          r_index(:) = int(r_scaled(:,i_part))  + 1

          ! DEBUG/CORRECT
          if(debug) then 
              if ( r_index(1) > nx ) then
                  print '(a,i6,a)',' index X out of bounds :',r_index(1),'Cycling ...' 
                  !            r_index(1) = nx
                  cycle
              end if
              if ( r_index(1) < 1 ) then
                  print '(a,i6,a)',' index X out of bounds :',r_index(1),'Cycling ...' 
                  !            r_index(1) =1 
                  cycle
              end if
              ! Z     
              if ( r_index(3) > nz ) then
                  print '(a,i6,a)',' index Z out of bounds :',r_index(3),'Cycling ...' 
                  !            r_index(1) = nx
                  cycle
              end if
              if ( r_index(3) < 1 ) then
                  print '(a,i6,a)',' index X out of bounds :',r_index(3),'Cycling ...' 
                  !            r_index(1) =1 
                  cycle
              end if
          end if
          !
          ! ------ Droplet    
          if(a_type(i_part)==3) then
              histo(r_index(1),1,r_index(3)) =  histo(r_index(1),1,r_index(3)) + 1. 
              histo_step(r_index(1),r_index(3)) =  histo_step(r_index(1),r_index(3)) + 1. 
          end if
          ! ------ Brush    
          !    if(a_type(i_part)==1.or.a_type(i_part)==2) then
          if(a_type(i_part) <  3 ) then  ! if brush
              histo_b(r_index(1),1,r_index(3)) =  histo_b(r_index(1),1,r_index(3)) + 1. 
              histo_step_b(r_index(1),r_index(3)) =  histo_step_b(r_index(1),r_index(3)) + 1. 
          end if

#ifdef FORCE_PROFILE
          if(a_type(i_part) == 1 ) then  ! if brush head 

              if(r_index(1)>nx.or.r_index(1)<1) print*,' BRUSH OUT OF BONDS !',i_part 

              histo_f(r_index(1)) =  histo_f(r_index(1)) + 1    

              f_index_x = (i_part-1)/n_mon + 1

              force_prof(r_index(1),:) =  force_prof(r_index(1),:) +  f_on_heads(:,f_index_x) 

          end if
#endif
      end do

! Squared mean values for density fluc.  profiles 
    
    do i = 1,nx
        do k = 1,nz
            histo2(i,1,k) = histo2(i,1,k) + histo_step(i,k)**2
            histo_b2(i,1,k) = histo_b2(i,1,k) + histo_step_b(i,k)**2
        enddo  !
    enddo

end select

end subroutine dens_prof
! 
!  ---- Velocity profile 2D for the droplet and brush 
!
   subroutine  vel_prof_2d (mode,histo,histo2,histo_b,histo_b2,n_box,r_box,dens_prof) ! 
     integer, intent(in) :: mode,n_box(3)
     real (kind=8), dimension(:,:,:) ,intent(inout) :: histo, histo_b, histo2,histo_b2
     real (kind=8), dimension(:,:,:) ,intent(in) :: dens_prof
     real (kind=8), intent(in) :: r_box(3)
     real (kind=8) ,allocatable,save :: r_scaled(:,:)
     integer,save :: i_time,nx,nz
     integer, allocatable, save :: icount(:,:) , icount_b(:,:)
     real (kind=8) , save :: dens_lim,inv_r_box(3)
     real (kind=8) :: max_vel,mod_vel,sig2
     integer :: i,j,k,i_part,r_index(3)
     logical , parameter :: debug=.false.

select case (mode)
case(0)    ! ---- Initialize 

    print "(/a/)", "  *   Initialising vel_prof_2d"
! Velocities     
    histo(:,:,:) = 0.0
    histo_b(:,:,:) = 0.0

    histo2(:,:,:) = 0.0
    histo_b2(:,:,:) = 0.0

! Counters

    allocate (icount(size(histo,dim=1),size(histo,dim=2) ) )
    allocate (icount_b(size(histo,dim=1),size(histo,dim=2)  ) )

    icount(:,:) = 0 
    icount_b(:,:) = 0 
    
    nx = size(histo,dim=1)
    nz = size(histo,dim=2)

    i_time = 0 
    print *, "nx,nz,n_part   = ",nx,nz,n_part
    
    allocate (r_scaled(3,n_part))

    r_scaled(1:3,:) = 0.0

    dens_lim = 0.5 ! print out only mean vels if the density is bigger than this param
    print '(/a,f5.2/)' ," *  Writing vel prof for density >= ", dens_lim

    inv_r_box(:) = 1./r_box(:)
    

    icount(:,:) = 0 
    icount_b(:,:) = 0 

case(1) ! --------- Normalization and writing

!here we write counters to file without normalizing    do i = 1, nx
!here we write counters to file without normalizing        do k= 1,nz
!here we write counters to file without normalizing
!here we write counters to file without normalizing        if ( icount(i,k) > 0
!) then
!here we write counters to file without normalizing          histo(i,k,:) =
!histo (i,k,:) /  real (icount(i,k)) ! the counters have the total number => no
!normalization with time
!here we write counters to file without normalizing          histo2(i,k,:) =
!histo2 (i,k,:)      /  real (icount(i,k))
!here we write counters to file without normalizing        end if
!here we write counters to file without normalizing        
!here we write counters to file without normalizing        if ( icount_b(i,k) >
!0 ) then
!here we write counters to file without normalizing          histo_b(i,k,:) =
!histo_b (i,k,:)  / real (icount_b(i,k))
!here we write counters to file without normalizing          histo_b2(i,k,:) =
!histo_b2 (i,k,:)  / real (icount_b(i,k))
!here we write counters to file without normalizing        end if
!here we write counters to file without normalizing
!here we write counters to file without normalizing         end do 
!here we write counters to file without normalizing    end do

!---  Write out droplet profile 

      print '(/a)', "DROP: Writing vel_prof_2d.mide ... " 
      

    open(unit=37,file="vel_prof_2d.mide",status="unknown")

! ----- Droplet velocity profile

    do k = 1,nz
           do i = 1,nx
!not here. Write all                if(dens_prof(i,1,k)> dens_lim) then ! only
!if the density is appreciable, write 

                write(37,'(6(f16.8,2x))')    dble(i-1)*r_box(1) + r_box(1)/2., &
                                             dble(k-1)*r_box(3) + r_box(3)/2., &
                                             histo(i,k,:),dble(icount(i,k))

!write all                                    else
!write all                 write(37,'(5f16.8)')    real(i-1)*r_box(1) +
!r_box(1)/2.,                        &
!write all                                         real(k-1)*r_box(3) +
!r_box(3)/2., 0.,0.,0.
!Write all                end if
           enddo
    enddo
    close(37)

!---  Write out brush profile 

      print '(a/)', "BRUSH:Writing vel_prof_2d_b.mide ... " 
    open(unit=37,file="vel_prof_2d_b.mide",status="unknown")
    
    !NOTE: The tot. vel. in each bin (without normalization ) and the value of
    !the counter are done.

! Header for automatic reading from tecplot
!not here         write(37,*) 'variables= "x","z","density"'
!not here         write(37,*) 'zone i= ',nz,' ,j= ',nx,' ,f=point'
!not here         write(37,*) 

! Brush 
    do i = 1,nx
                do k = 1,nz
                write(37,'(6(f16.8,x))')    dble(i-1)*r_box(1) + r_box(1)/2.,  &
                                        dble(k-1)*r_box(3) + r_box(3)/2., histo_b(i,k,:),dble(icount_b(i,k))

                enddo
    enddo
    close(37)

! --- Sigma^2(v) profile 


      print '(/a)', "DROP: Writing sig2_vel_2d.mide ... " 
      

    open(unit=37,file="sig2_vel_2d.mide",status="unknown")

! Header for automatic reading from tecplot

       write(37,*) 'variables= "x","z","dens"'
       write(37,*) 'zone i= ',nz,' ,j= ',nx,' ,f=point'
       write(37,*) 

! droplet velocity profile
        do k = 1,nz
                    do i = 1,nx
                    if(dens_prof(i,1,k) > dens_lim) then ! only if the density is appreciable, write 
                        
                        sig2= sum ( histo2(i,k,:) - histo(i,k,:)**2 )
                        
                    write(37,'(3e16.8)')    dble(i-1)*r_box(1) + r_box(1)/2., &
                                            dble(k-1)*r_box(3) + r_box(3)/2., sig2
                                        else
                    write(37,'(3e16.8)')    dble(i-1)*r_box(1) + r_box(1)/2., &
                                            dble(k-1)*r_box(3) + r_box(3)/2., 0.
                    end if
                    enddo
        enddo
        close(37)

      case(3) ! ------  Order N algorithm  ! WORKS ! 

! time counter
      i_time = i_time + 1
      
! change scale so that int() means find the right binning box for the particle      

          r_scaled(1,:) = r0_unfold(1,:)*inv_r_box(1)
          r_scaled(2,:) = r0_unfold(2,:)*inv_r_box(2)
          r_scaled(3,:) = r0_unfold(3,:)*inv_r_box(3)
  
!    ------ Brush    
    do i_part = 1, part_init_d      

          r_index(:) = int(r_scaled(:,i_part))  + 1

          if(debug) then
             if (r_index(1) > size(histo,dim=1) ) print *, " Brush: Bounds exceeded !  ",r_index(1)
             if (r_index(3) > size(histo,dim=2) ) print *, " Brush: Bounds exceeded !  ",r_index(3)
          end if


          histo_b(r_index(1),r_index(3),:) =  histo_b(r_index(1),r_index(3),:) + v(:,i_part)   
          histo_b2(r_index(1),r_index(3),:) =  histo_b2(r_index(1),r_index(3),:) + v(:,i_part)**2   
!    counter
          icount_b(r_index(1),r_index(3)) = icount_b(r_index(1),r_index(3)) + 1

      end do

!    ------- Drops

    do i_part = part_init_d+1,n_part

          r_index(:) = int(r_scaled(:,i_part))  + 1

          if(debug) then
             if (r_index(1) > size(histo,dim=1) ) print *, " Drops: Bounds exceeded !  ",r_index(1)
             if (r_index(3) > size(histo,dim=2) ) print *, " Drops: Bounds exceeded !  ",r_index(3)
          end if

          histo(r_index(1),r_index(3),:) =  histo(r_index(1),r_index(3),:) + v(:,i_part)   
          histo2(r_index(1),r_index(3),:) =  histo2(r_index(1),r_index(3),:) + v(:,i_part)**2   
!    counter

          icount(r_index(1),r_index(3)) = icount(r_index(1),r_index(3)) + 1
    end do


end select

   end subroutine  vel_prof_2d 

! X coordinate refolding routine: paste the droplet in x dir
! WARN: assumes that all the free chains have the same length

subroutine coor_refold(mode)
       integer, intent(in) :: mode
       integer :: i
       
       select case(mode)
       case(1)
      ! Warning: twice the box  
           do i=1,n_part
           if (r0_unfold(1,i) > 2*boundary(1) ) r0_unfold(1,i) = r0_unfold(1,i) - boundary(1) 
           if (r0_unfold(1,i) < 0.  ) r0_unfold(1,i) = r0_unfold(1,i) + boundary(1) 
           end do 
       case(2)
       
           do i=1,n_part
           if (r0_unfold(2,i) > boundary(2) ) r0_unfold(2,i) = r0_unfold(2,i) - boundary(2) 
           if (r0_unfold(2,i) < 0.  ) r0_unfold(2,i) = r0_unfold(2,i) + boundary(2) 
           end do 
       case(3)
       
           do i=1,n_part
           if (r0_unfold(3,i) > boundary(3) ) r0_unfold(3,i) = r0_unfold(3,i) - boundary(3) 
           if (r0_unfold(3,i) < 0.  ) r0_unfold(3,i) = r0_unfold(3,i) + boundary(3) 
           end do 

       case(4)
      ! Warning: Lx = boundary(1)  
           do i=1,n_part
               if (r0_unfold(1,i) > boundary(1) ) r0_unfold(1,i) = r0_unfold(1,i) - boundary(1) 
               if (r0_unfold(1,i) < 0.  ) r0_unfold(1,i) = r0_unfold(1,i) + boundary(1) 
           end do 

       
       end select

end subroutine coor_refold

#   endif /*SYSTEM==1; droplets*/


!
!  ENDS MODULE        
   end module functions
