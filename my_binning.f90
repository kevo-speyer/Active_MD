!   Copied from Cem Servantie implemmentation
!  - Does not use newton 3 law. Goes over all the particles 
!  - taken from Plimpton's implementations [ref: J. Comp. Phys 117, 1, '95]

    subroutine my_binning()
#include 'control_simulation.h'     
    use commons 
#ifdef MPI
    use mpi_commons
#endif
    implicit none
    logical,parameter :: binning_debug=.false.
    
!    integer :: binpnt(n_max_cell),bin(n_max_cell)
    integer :: ix,iy,iz,ib
    integer :: ixx,iyy,izz,k
    integer :: n1,n2
    real (kind=8) :: r0_tmp(3),dr(3),dr_2
    
! set to zero 
    
    ff_list(:,1:n_part) = 0 
    binpnt(:) = 0
    bin(:) = 0
    
    do i_part=1,n_part !n_mon_tot
    
        ix = int ( r0(1,i_part) / r_bin_x  )  
        iy = int ( r0(2,i_part) / r_bin_y  ) 
        iz = int ( r0(3,i_part) / r_bin_z  ) 
    
        ib = iz*n_bin_x*n_bin_y + iy*n_bin_x+ ix + 1   
    
        bin(i_part) = binpnt(ib)  
        binpnt(ib) = i_part
! check it
! CEM: not sure if doing it always
        if(binning_debug) then
            if (ib<1 .or. ib>n_cell) then
                print '(a,4i6)', "ERROR in BINNING: i_part,ix,iy,iz=", i_part,ix,iy,iz
                print '(a,3f10.5)', "ERROR in BINNING: part. pos= ", r0(:,i_part)
            end if
        end if
    end do
    
!PARALLEL: here calls to the interval routine and gets n_min n_max for the
!current processor
    
! CAN BE PARALLEL
#ifdef MPI
        call interval(my_rank,n1,n2)         
!        print *,"N1 N2=",n1,n2
        do i_part= n1,n2
#else 
        do i_part=1,n_part
#endif        

    
        r0_tmp(:) = r0(:,i_part)
        i_type = a_type (i_part) ! identify i_part particle
    
        ixx = int ( r0(1,i_part) / r_bin_x )  
        iyy = int ( r0(2,i_part) / r_bin_y ) 
        izz = int ( r0(3,i_part) / r_bin_z ) 
    
! go over all the neighbouring cells (27) without Newton 3
    
        do k = 0,26
    
            ix = ixx + mod(k,3) - 1 
            if (ix <  0) ix = n_bin_x -1
            if (ix == n_bin_x) ix = 0
    
            iy = iyy + mod(k/3,3) - 1 
            if (iy <  0) iy = n_bin_y -1
            if (iy == n_bin_y) iy = 0
    
            iz = izz + k/9 - 1 
            if (iz <  0) iz = n_bin_z -1
            if (iz == n_bin_z) iz = 0
    
            ib = iz*n_bin_x*n_bin_y + iy*n_bin_x+ ix + 1   
    
            j_part = binpnt(ib)
    
10          if(j_part /= 0) then
    
                j_type = a_type(j_part)
                dr(:) = r0_tmp(:) - r0(:,j_part)
    
! Apply PBC (x and y)
#       if SYMMETRY == 0 
                do i_dim = 1,n_dim-1
                    dr(i_dim) = dr(i_dim) - boundary(i_dim)*int(2*dr(i_dim)*inv_boundary(i_dim))
                end do
#       elif SYMMETRY == 1
                do i_dim = 1,n_dim
                    dr(i_dim) = dr(i_dim) - boundary(i_dim)*int(2*dr(i_dim)*inv_boundary(i_dim))
                end do
#       endif
    
                dr_2 =  dr(1)**2+dr(2)**2+dr(3)**2 
    
! if it is inside the radius, add up to the list
    
                if ( dr_2 < range_2(i_type,j_type)+skin_2 .and. i_part /= j_part  ) then
                    ff_list(0,i_part) =  ff_list(0,i_part) + 1
                    ff_list(ff_list(0,i_part),i_part ) = j_part
                end if
                j_part = bin (j_part)
                goto 10 
            end if  
        end do  ! loop on neighboring cells 

        if(binning_debug) then
            if(ff_list(0,i_part) > n_bin_fl) then
                print '(/a,3i5)'," ERROR: Neighbor list too big: current,max,i_part= ", & 
                    ff_list(0,i_part),n_bin_fl,i_part
                print '(/a,1i5,3f12.8/)',"i_part,r0(i_part)= ", i_part,r0(:,i_part)
             !   print *,n_bin_fl
                stop
            end if
        end if
    
    end do  ! i_part
    
! Update check_skin variables (inherited from binning.f90)
    
    f_skin = 0
    r0_old(:,1:n_part) = r0(:,1:n_part)
    mic_old(:,1:n_part) = mic_count(:,1:n_part)
    
    
    if(binning_debug) then 
!        print '(/a/)',"DEBUG ON in my_binning.f90" 
!        print *, ff_list(1:3,:) ; stop
    end if
    
    
    end subroutine my_binning
