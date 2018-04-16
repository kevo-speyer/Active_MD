!
!-------- Checks if the binning will be done in this step 
!
      subroutine check_skin
      use commons ; implicit none
      logical, parameter :: debug=.false.

              do i_part = 1,n_part
                  do i_dim = 1,n_dim
                      delta_r(i_dim) =                                                     &
                          r0_old(i_dim,i_part)+mic_old(i_dim,i_part)*boundary(i_dim)       &
                          -r0(i_dim,i_part)-mic_count(i_dim,i_part)*boundary(i_dim)
                  end do
                  r_dummy = 0.
                  do i_dim = 1,n_dim
                      r_dummy = r_dummy + delta_r(i_dim)**2
                  end do
                  if (r_dummy.gt.(skin/2.)**2) then
                      f_skin = 1
                  end if
              end do
              if(debug .and. f_skin == 1) print '(a,i6)',"Doing binning in step",i_time  
!
      end subroutine check_skin
