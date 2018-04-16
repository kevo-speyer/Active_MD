! Put debugging entrances with parameter binning debug
      subroutine binning
#include 'control_simulation.h'
      use commons ; implicit none
      logical,parameter :: binning_debug=.true.

!cla: I think that is the linked-list algorithm (Smit & Frenkel, 1996, p. 368 )
!---  zero the binning entrances

    do i_bin_z = 0,n_bin_z-1
        do i_bin_y = 0,n_bin_y-1
            do i_bin_x = 0,n_bin_x-1
                if(binning_debug) then
                    n_fcell_max = max(bin_fluid(i_bin_x,i_bin_y,i_bin_z,0),n_fcell_max)
                end if
                bin_fluid(i_bin_x,i_bin_y,i_bin_z,0:n_bin_fl) = 0 
            end do
        end do
    end do


!---- bin the particles into the boxes

!ori    do i_part = 1,n_mon_tot
    do i_part = 1,n_part

        i_bin_x = mod(int(r0(1,i_part)*inv_r_bin_x),n_bin_x)
        i_bin_y = mod(int(r0(2,i_part)*inv_r_bin_y),n_bin_y)
        i_bin_z = mod(int(r0(3,i_part)*inv_r_bin_z),n_bin_z)
!debug
        if(binning_debug) then
            if (i_bin_x>n_bin_x-1.or.i_bin_y>n_bin_y-1.or.i_bin_z>n_bin_z-1) then
                print'(a,3i6)',"WARN: some i_bin_? bigger that n_bin_? limits exceeded: i_bin_[xyz]", &
                i_bin_x,i_bin_y,i_bin_z
            end if
        end if 

        bin_fluid(i_bin_x,i_bin_y,i_bin_z,0) = bin_fluid(i_bin_x,i_bin_y,i_bin_z,0) + 1

        if(binning_debug) then 
!            print *,bin_fluid(i_bin_x,i_bin_y,i_bin_z,0)
            if ( bin_fluid(i_bin_x,i_bin_y,i_bin_z,0) > n_bin_fl+ 1 ) then
                print '(a,4i6)', "WARN: max dim in bin_fluid excedeed: i_bin_[xyz], bin_fluid", & 
                i_bin_x,i_bin_y,i_bin_z,bin_fluid(i_bin_x,i_bin_y,i_bin_z,0) 
            end if
        end if

        bin_fluid(i_bin_x,i_bin_y,i_bin_z,bin_fluid(i_bin_x,i_bin_y,i_bin_z,0)) = i_part

    end do

!
!---- generate lists
!
    !old do i_part = 1,n_mon_tot
     do i_part = 1,n_part
        n_neigh_fl_max = max(n_neigh_fl_max,ff_list(0,i_part))
        ff_list(0,i_part) = 0
    end do
    
!
!---- Build up fluid fluid list and DPD list. They can differ now.
!
    do i_bin_z = 0,n_bin_z-1
        do i_bin_y = 0,n_bin_y-1
            do i_bin_x = 0,n_bin_x-1
                n1_bin_part = bin_fluid(i_bin_x,i_bin_y,i_bin_z,0)
                do i1_bin_part = 1,n1_bin_part
                    i_part = bin_fluid(i_bin_x,i_bin_y,i_bin_z,i1_bin_part)
                    i_type =a_type(i_part)
                    do j_bin_x = i_bin_x-1,i_bin_x+delta_bin_x
                        k_bin_x = mod(j_bin_x+n_bin_x,n_bin_x)
                        do j_bin_y = i_bin_y-1,i_bin_y+delta_bin_y
                            k_bin_y = mod(j_bin_y+n_bin_y,n_bin_y)
                            do j_bin_z = i_bin_z-1,i_bin_z+delta_bin_z
                                k_bin_z = mod(j_bin_z+n_bin_z,n_bin_z)
                                n2_bin_part = bin_fluid(k_bin_x,k_bin_y,k_bin_z,0)
                                do i2_bin_part = 1,n2_bin_part
                                    j_part = bin_fluid(k_bin_x,k_bin_y,k_bin_z,i2_bin_part)

                                    if(i_part.lt.j_part) then
                                        j_type =a_type(j_part)
!HPC (?)
                                        delta_r(1) = r0(1,i_part) - r0(1,j_part)
                                        delta_r(2) = r0(2,i_part) - r0(2,j_part)
                                        delta_r(3) = r0(3,i_part) - r0(3,j_part)

!-----  PBC 
                                        delta_r(1) = delta_r(1) - boundary(1)*int(2*delta_r(1)*inv_boundary(1))
                                        delta_r(2) = delta_r(2) - boundary(2)*int(2*delta_r(2)*inv_boundary(2))
#if SYMMETRY == 1
                                        delta_r(3) = delta_r(3) - boundary(3)*int(2*delta_r(3)*inv_boundary(3))
#endif
                                        r_2 = 0.
! HPC
                                        r_2 = r_2 + delta_r(1)*delta_r(1)
                                        r_2 = r_2 + delta_r(2)*delta_r(2)
                                        r_2 = r_2 + delta_r(3)*delta_r(3)

!-----  Check if distance is inside the cutoff sphere
#ifdef DPD_CUT_OFF
                                        if(r_2.lt.r_cut_max_2+skin_2) then
#else
                                        if(r_2.lt.range_2(i_type,j_type)+skin_2) then
#endif                                            
                                            if(ff_list(0,i_part).le.ff_list(0,j_part)) then
                                                ff_list(0,i_part) = ff_list(0,i_part) + 1
                                                ff_list(ff_list(0,i_part),i_part) = j_part
                                            else
                                                ff_list(0,j_part) = ff_list(0,j_part) + 1
                                                ff_list(ff_list(0,j_part),j_part) = i_part
                                            end if
                                        end if
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do !ends fluid-fluid binning


      f_skin = 0
      r0_old(:,1:n_part) = r0(:,1:n_part)
      mic_old(:,1:n_part) = mic_count(:,1:n_part)

!if(binning_debug) print *,"Maximum neighbors:", maxval(bin_fluid(:,:,:,0))
!pause
end subroutine binning
