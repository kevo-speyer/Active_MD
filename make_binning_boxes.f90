! This routine sets up the size and number of binning boxes
! using the cutoff distance of the simulated system
! This will be used to allocate vectors for del linked-cell list algorihtm
! implemmented in the "binning" routine.
subroutine make_binning_boxes()
    use commons
#   include 'control_simulation.h'    
implicit none
logical, parameter :: debug = .true.

!
!---- Initialize binning
!
      i_dummy = n_bin_x
      if(i_dummy.ge.3) then
      delta_bin_x = 1
      else if (i_dummy.eq.2) then
      delta_bin_x = 0
      else if (i_dummy.eq.1) then
      delta_bin_x = -1
      end if
      i_dummy = n_bin_y
      if(i_dummy.ge.3) then
          delta_bin_y = 1
         else if (i_dummy.eq.2) then
          delta_bin_y = 0
         else if (i_dummy.eq.1) then
          delta_bin_y = -1
      end if
      i_dummy = n_bin_z
      if(i_dummy.ge.3) then
          delta_bin_z = 1
         else if (i_dummy.eq.2) then
          delta_bin_z = 0
         else if (i_dummy.eq.1) then
          delta_bin_z = -1
      end if
!
! --- Build up binning boxes
!
binning_box_skin  = 0.4   

n_bin_x = int( boundary(1)/(r_cut_max + binning_box_skin)  )   
n_bin_y = int( boundary(2)/(r_cut_max + binning_box_skin)  )   
n_bin_z = int( boundary(3)/(r_cut_max + binning_box_skin)  )   

r_bin_x = boundary(1)/dble(n_bin_x)      !ori: r_bin_x = boundary(1)/n_bin_x     
r_bin_y = boundary(2)/dble(n_bin_y)      !ori: r_bin_y = boundary(2)/n_bin_y     
r_bin_z = boundary(3)/dble(n_bin_z)      !ori: r_bin_z = 2.5 !boundary(3)/n_bin_z

! HPC
inv_r_bin_x =   1./r_bin_x
inv_r_bin_y =   1./r_bin_y
inv_r_bin_z =   1./r_bin_z
!write(*,*) boundary(:), r_cut_max, binning_box_skin

!
! Allocation of binning variables
!
allocate ( bin_twall(0:n_bin_x-1,0:n_bin_y-1,0:n_bin_z-1,0:n_bin_wa),   &
           bin_bwall(0:n_bin_x-1,0:n_bin_y-1,0:n_bin_z-1,0:n_bin_wa),   &
           bin_fluid(0:n_bin_x-1,0:n_bin_y-1,0:n_bin_z-1,0:n_bin_fl)  )

! Allocation of variables for my_binning 

    n_cell = n_bin_x*n_bin_y*n_bin_z 

    allocate(binpnt(n_bin_x*n_bin_y*n_bin_z+n_bin_x*n_bin_y+n_bin_x+1),bin(n_part))

!NCELLX*NCELLY*NCELLZ+NCELLY*NCELLX+NCELLX+1
print '(/a,3i6)',     "   - Number of binning box in each coor. = ",n_bin_x,n_bin_y,n_bin_z
print '(a,3f7.3/)',    "   - Binning box dimensions =              ",r_bin_x,r_bin_y,r_bin_z


!
! Write out if binning is not coherent with binning boxes
!
        write(20,*)
        write(20,'(3i15,a)') n_bin_x,n_bin_y,n_bin_z,"  binning boxes"
        write(20,'(3f15.6,a)') r_bin_x,r_bin_y,r_bin_z,"  binning distances"
        do i_type = 1,n_type
            do j_type = i_type,n_type
                if((range_2(j_type,i_type)+skin_2.gt.r_bin_x**2).or.             &
                    (range_2(j_type,i_type)+skin_2.gt.r_bin_y**2).or.            &
                    (range_2(j_type,i_type)+skin_2.gt.r_bin_z**2)) then
                write(*,*) j_type,i_type," binning problematic"
                write(*,'(3f15.7)') sqrt(range_2(j_type,i_type)),r_bin_x,r_bin_y,r_bin_z
            end if
        end do
    end do


end subroutine make_binning_boxes
