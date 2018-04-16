      subroutine wall_wall(inter_type)
      use commons ; implicit none
      integer, intent(in) :: inter_type   
      v_wall_wall = 0.
!
select case (inter_type)
! now is done in fluid fluid 
case(1)         ! LJ wall particles interaction
!obs    i_type = n_type
!obs    j_type = n_type
!obs    do i_wall = 1,n_wall/2
!obs        i_part = i_wall + n_mon_tot
!obs        i_dummy = ww_list(0,i_wall)
!obs        do i_neigh = 1,i_dummy
!obs            j_part = ww_list(i_neigh,i_wall)
!obs            do i_dim = 1,n_dim
!obs                delta_r(i_dim) = r0(i_dim,i_part) - r0(i_dim,j_part)
!obs            end do
!obs            !-----  get boundaries right
!obs            do i_dim = 1,n_dim-1
!obs                delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*             &
!obs                    &   int(2*delta_r(i_dim)/boundary(i_dim))
!obs            end do
!obs            r_2 = 0.
!obs            do i_dim = 1,n_dim
!obs                r_2 = r_2 + delta_r(i_dim)**2
!obs            end do
!obs            !-----  check whether interaction takes place
!obs            if(r_2.lt.range_2(i_type,j_type)) then
!obs                r_2 = max(r_2,r_2_min)!(i_type, j_type))
!obs                r_6 = (sigma_2(i_type,j_type)/r_2)**3
!obs                r_12 = r_6**2
!obs                pot_loc = (r_12-r_6) - e_shift(i_type,j_type)
!obs                v_wall_wall = v_wall_wall + epsil(i_type,j_type)*pot_loc
!obs                r_dummy = epsil(i_type,j_type)*(-12*r_12+6*r_6)/r_2
!obs                do i_dim = 1,n_dim
!obs                    force_loc(i_dim) = r_dummy*delta_r(i_dim)
!obs                    force(i_dim,i_part) = force(i_dim,i_part) - force_loc(i_dim)
!obs                    force(i_dim,j_part) = force(i_dim,j_part) + force_loc(i_dim)
!obs                end do
!obs            end if
!obs        end do
!obs    end do
! only interaction between two walls ! IS IT OK ??
      case (2) 
#       if WALL==2      
          v_wall_wall = 2.0*( abs(a_w)*(sigma_w/z_space_wall)**9 - a_w*(sigma_w/z_space_wall)**3 ) 
#       endif    
      case(3) ! Only bottom wall has 9-3 potential and top wall is a hard wall 
#       if WALL==3      
          v_wall_wall = 1.0*( abs(a_w)*(sigma_w/z_space_wall)**9 - a_w*(sigma_w/z_space_wall)**3 ) 
#       endif    
      case(4) ! Top and bottom walls are both hard walls
          ! NO force between walls  
          continue



      case default
          print*," inter_type must be 1: LJ-wall or 2:no interaction in walls, 3 or 4 "
          print*,"Change it ! Stopping here."
          stop
      end select

end subroutine wall_wall
