!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
! HOW TO USE THE LINKED CELL LIST:                                  !
! === == === === ====== ==== ====                                   !
!                                                                   !
! (  LET'S GET READY TO RESPA EDITION    )                          !
!                                                                   !
! 1) call the routine "count_cells", with the proper arguments.     !
!    This routine should only be called once.                       !
!                                                                   !
! 2) call the routine "neigh_list", to generate linked lists:       !
!    part_in_cell, lpart_in_cell, r_nei, l_nei.                     !
!    This routine should only be called once.                       !
!                                                                   !
! 3) call get_cell_neigh_list to generate lists of all neighboring  !
!    cells for each cell.                                           !
!    This routine should only be called once.                       !
!                                                                   !
! 4) to loop particles in a cell:                                   !
!   i_part = part_in_cell(i_cell)                                   !
!   do while(i_part .ne. 0)                                         !
!       !Do whatever with r0(:,i_part)                              !
!       i_part = r_nei(i_part) !gets next particle in this cell     !
!   end do                                                          !
!                                                                   !
! 5) To loop between possible interacting particles of              !
!       a particle (i_part) :                                       !
!   a) call get_cell, to obtain the cell where i_part               !
!        is located (i_cell)                                        !
!   b) loop between neighboring cells of i_cell neig_cells(i_cell,:)!
!   c) loop over all particles in those cells (see 4), above)       !
!   do j=1,3**n_dim-1                                               !
!       j_cell = neig_cells(i_cell,j)                               !
!       j_part = part_in_cell(j_cell)                               !
!       do while(j_part .ne. 0)                                     ! 
!           !Do whatever with r0(:,j_part) and r0(:,i_part)         !  
!           j_part = r_nei(j_part) !gets next particle in this cell !    
!       end do                                                      !
!   end do                                                          !
!   d) loop over particles in the same cell also (i_cell)           !
!   j_part = part_in_cell(i_cell)                                   !
!   do while(j_part .ne. 0)                                         !
!       if(j_part .eq. i_part) then                                 !
!           j_part = r_nei(i_part)                                  !
!           cycle    ! don't count auto-interactions                !
!       end if                                                      !
!       !Do whatever with r0(:,i_part) and r0(:,j_part)             !
!       j_part = r_nei(j_part) !gets next particle in this cell     !
!   end do                                                          !
!                                                                   !
! 6) Check if particle i_part leaves a cell and update list:        !
!    a) Before updating position, get it's cell (old_cell)          !
!    b) call subroutine  "update_part_cell"                         !
!    This will automatically check if particle changes cell,        !
!    and will update linked lists  and cell list                    !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine neigh_list(r0, a_type, n_part, n_dim, n_cells, n_cells_tot, boundary, inv_l_cell, part_in_cell, lpart_in_cell, r_nei, l_nei, cell_of_part)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to make cell lists in a linked way                        !
! part_in_cell(i_cell) gives the index of the first                 !
! particle in the cell.                                             !
! lpart_in_cell(i_cell) gives the index of the                      !
! last particle in the cell                                         !
! r_nei(i_part) gives the next particle in the same cell.           !
! If there are no more particles in this cell,                      !
! then r_nei(i_part) = 0                                            !
! l_nei(j_part) gives the previous particle in the cell. If this is !
! the first particle in the cell, then l_nei(j_part) = 0            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(kind=8), intent(in) :: inv_l_cell(n_dim), r0(n_dim,n_part), boundary(n_dim)
integer, intent(in) :: n_part, n_dim,n_cells(n_dim), n_cells_tot, a_type(n_part)
integer, intent(out) :: part_in_cell(0:n_cells_tot), lpart_in_cell(0:n_cells_tot), r_nei(n_part), l_nei(n_part), cell_of_part(n_part)
integer :: i_part, j_part, p_cell, i_dim

!!!!!!!!! DO LINKED LIST FROM SCRATCH !!!!!!!!!!!!!!!!
part_in_cell(0:n_cells_tot) = 0
lpart_in_cell(0:n_cells_tot) = 0
l_nei(:) = 0
r_nei(:) = 0
do i_part = 1, n_part

!! First check that the  particle is in the simulation box
    !print*,"check_bound for i_part",i_part
call check_bound(r0(:,i_part), n_dim, boundary)

!print*,"Doing the linked list for particle ",  i_part
    !First get the cell index of the particle. Cell index is between 1 and
    !n_cells_tot
    call get_cell(r0(:,i_part), n_dim, n_cells, inv_l_cell, p_cell)
    cell_of_part(i_part) = p_cell
    !if(p_cell.gt.n_cells_tot) print*, "ERROR cell out of bounds",r0(:,i_part),p_cell 
    ! If Cell is empty, link cell with particle i_part
    if( part_in_cell(p_cell) .eq. 0 ) then
        part_in_cell(p_cell) = i_part
        lpart_in_cell(p_cell) = i_part
    else !If cell has already a particle in 
        if ( a_type(i_part) .ne. 3) then ! insert particl by the left (head or first part.)
            j_part = part_in_cell(p_cell)
            part_in_cell(p_cell) = i_part
            r_nei(i_part) = j_part
            l_nei(i_part) = 0 !Done already
            l_nei(j_part) = i_part
        else                        ! insert particle by the right (tail or last part.)
            j_part = lpart_in_cell(p_cell)
            lpart_in_cell(p_cell) = i_part
            l_nei(i_part) = j_part
            r_nei(i_part) = 0 !Done already
            r_nei(j_part) = i_part
        end if
    end if
end do

! To loop faster, a phantom cell (cell_index=0) is created. This cell never has
! any particles. This is usefull if Periodic Boundary Conditions are NOT applied
! in z to loop over cells that are in the cell list, but are not physically next
! to each other. (See subroutine make_cell_list2 to see implementatio to see
! implementationn)
part_in_cell(0) = 0
lpart_in_cell(0) = 0

end subroutine

subroutine check_bound(r0_part, n_dim, boundary)
implicit none
real(kind=8), intent(in) :: r0_part(n_dim), boundary(n_dim)
integer, intent(in) :: n_dim !, n_part
integer :: i_dm, j_dm

!check if in particle are inside the simulation box

!do i_dm = 1; n_part
    do j_dm = 1, n_dim
        if( (r0_part(j_dm).lt.0.) .or. (r0_part(j_dm).gt.boundary(j_dm)) ) then
            print*, "ERROR in conf_old. Particle outside boundaries"
            print*, "r0(:,i_part) = ",r0_part(:)
            print*, "boundaies = ", boundary(:) 
            stop
        end if
    end do
!end do

end subroutine check_bound

subroutine make_binning(bin_case, n_dim, boundary, r_cut, l_cell, inv_l_cell, n_cells, n_cells_tot, n_nei)
implicit none
integer, intent(in) :: n_dim, bin_case
real(kind=8), intent(in) :: boundary(n_dim), r_cut
real(kind=8), intent(out) :: l_cell(n_dim), inv_l_cell(n_dim)
integer, intent(out) :: n_cells(n_dim), n_cells_tot, n_nei
integer :: i_dim

select case (bin_case)
case(1) !l_cell = r_cut ! 13 neighbors per cell 
n_cells(:) = int( boundary(:) / r_cut ) !count number of cells in each direction

n_nei = int( ( 3**n_dim - 1 ) / 2 )

case(2) !l_cell = r_cut / 2 ! 62 neighbors per cell 
n_cells(:) = int( boundary(:) / (r_cut / 2.) )

n_nei = int( ( 5**n_dim - 1 ) / 2 ) ! if l_cell = r_cut / 2

end select

l_cell(:) =  boundary(:) / float(n_cells(:)) ! take a different value of l_cell in each
inv_l_cell = 1 / l_cell
!Set  total number of  cells in the system
n_cells_tot = 1
do i_dim = 1, n_dim
    n_cells_tot = n_cells_tot * n_cells(i_dim)
end do

end subroutine make_binning

!subroutine count_cells(n_dim, boundary, l_cell, n_cells, n_cells_tot)
!implicit none
!!OBSOLETE NOW.
!! USE make_binning istead
!!This routine does the binning, and output is number of cells in each direction
!!and total number of cells
!integer, intent(in) :: n_dim ! dimension of space 
!real(kind=8), intent(in) :: boundary(n_dim), l_cell ! length of simulation box
!                                        !in each direction, linear length of a
!                                        !cell
!integer, intent(out) :: n_cells(n_dim), n_cells_tot !number of cells in each
!                                        !direction, total number fo cells
!integer :: i_dim
!!Cells are cubes of length l_cell
!! l_cell should be equal to r_cut
!!Set number of cells in each direction 
!n_cells(:) = int( boundary(:) / l_cell ) + 1
!
!!Set  total number of  cells in the system
!n_cells_tot = 1
!do i_dim = 1, n_dim
!    n_cells_tot = n_cells_tot * n_cells(i_dim)
!end do
!
!end subroutine

subroutine get_cell(r0_part, n_dim, m_cells, inv_l_cell, p_cell)
!Get the cell number of particle r0_part(:)
implicit none
real(kind=8), intent(in) :: r0_part(n_dim), inv_l_cell(n_dim)
integer, intent(in) :: n_dim, m_cells(n_dim)
integer, intent(out) :: p_cell
integer :: pf, i_dim
!debug
!integer :: n_cells_tot
!n_cells_tot = m_cells(1) * m_cells(2) * m_cells(3)

pf = 1
p_cell = 1
    do i_dim = 1, n_dim
        p_cell = p_cell + pf * int( r0_part(i_dim) * inv_l_cell(i_dim) )
        pf = pf * m_cells(i_dim) 
    end do

    !debug
    !if ((p_cell.gt.n_cells_tot).or.(p_cell.lt.1)) then
    !    print*, "Error particle out of cell"
    !    print*,"cell = ,",p_cell,"n_cells_tot",n_cells_tot
    !    print*, "r0 = ",r0_part(:)
    !    print*, "n_cells", m_cells(:)
    !    print*, "l_cell", l_cell
    !end if

end subroutine

subroutine update_part_cell(r0_part, i_part, n_dim, n_cells, n_cells_tot, inv_l_cell, n_part, part_in_cell, lpart_in_cell, r_nei, l_nei, h_or_t, cell_of_part)
!This routine updates the linked lists of the cells, if a particle leaves a cell

!! IN VARIABLES
!r0_part  is the coordinates of the position of a particle
! i_part  is the index of the particle
! n_dim   is the number of space dimnesions of the system
! n_cells is an array which contains the number of cells in each of the n_dim directions
! n_cells_tot is the total number of cells in the system
! l_cell  is the length of the cell box
! old_cell is the cell where the particle (i_part) was before updating positions
! h_or_t defines if the particle is added as the first particle or last particle of the linked list of that cell
! cell_of_part is an array which gives the cell in which each particle is located

!! IN OUT VARIABLES
! part_in_cell is an array that stores the index of the first particle for each cell or 0 if the cell is empty
! lpart_in_cell is an array that stores the index of the last particle for each cell or 0 if the cell is empty
! r_nei is the right cell neighbour of each particle
! l_nei is the left cell neighbour of each particle
! cell_of_part is the array that give the corresponding cell for each particle

implicit none
real(kind=8), intent(in) :: r0_part(n_dim), inv_l_cell(n_dim)
integer, intent(in) :: n_dim, n_cells(n_dim), i_part, n_cells_tot, n_part, h_or_t
integer, intent(inout) :: part_in_cell(0:n_cells_tot), lpart_in_cell(0:n_cells_tot), r_nei(n_part), l_nei(n_part), cell_of_part(n_part)
integer ::  new_cell, j_part, old_cell

old_cell = cell_of_part(i_part)

!get new cell
call get_cell(r0_part, n_dim, n_cells, inv_l_cell, new_cell)

if( new_cell .ne. old_cell ) then !Update cell list
    cell_of_part(i_part) = new_cell

    !Remove particle from old cell.
    if( lpart_in_cell(old_cell) .eq. i_part ) then ! if it's the las particle in the cell, relink lpart_in_cell to
        lpart_in_cell(old_cell) = l_nei( i_part )   ! an the left particle of i_part
    end if

    if( part_in_cell(old_cell) .eq. i_part ) then ! if it's the first particle in the cell, relink part_in_cell to
       part_in_cell(old_cell) = r_nei(i_part)       ! the particle to the right of i_part
    end if

! Re-link the nighbors that stay in this cell
if( l_nei(i_part) .ne. 0 ) r_nei( l_nei(i_part) ) = r_nei( i_part )

if( r_nei(i_part) .ne. 0 ) l_nei( r_nei(i_part) ) = l_nei( i_part ) 
   
   
    !Add particle to new cell 
    if( part_in_cell(new_cell) .eq. 0 ) then !If there are no particles in this cell
        part_in_cell(new_cell) = i_part !Put particle in new cell
        lpart_in_cell(new_cell) = i_part
        l_nei(i_part) = 0
        r_nei(i_part) = 0
    else !If cell has already a particle in
        if( h_or_t .eq. 1 ) then ! if the particle should be added by the left (head or first)
            j_part = part_in_cell(new_cell)
            part_in_cell(new_cell) = i_part
            r_nei(i_part) = j_part
            l_nei(i_part) = 0 !Done already
            l_nei(j_part) = i_part
        else if( h_or_t .eq. 2 ) then ! if the particle should be added by the right (tail or last)
            j_part = lpart_in_cell(new_cell)
            lpart_in_cell(new_cell) = i_part
            r_nei(i_part) = 0
            l_nei(i_part) = j_part !Done already
            r_nei(j_part) = i_part           
        end if
    end if
end if

end subroutine

subroutine get_cell_neigh_list(n_cells_tot, n_cells, n_dim, n_nei,n_nei_tot, cell_neigh_ls, cell_neigh_ls_tot)
!Generate a list of neighboring cells for each cell
implicit none
integer, intent(in) :: n_dim, n_cells_tot, n_cells(n_dim), n_nei
integer, intent(out) :: cell_neigh_ls(n_cells_tot,n_nei), cell_neigh_ls_tot(n_cells_tot,n_nei_tot)
integer :: i_cell, n_nei_tot

do i_cell = 1, n_cells_tot
    call make_cell_list2(i_cell, n_cells, n_dim, n_nei, cell_neigh_ls(i_cell,:))
    ! Only half the neighbors are given, in order to loop over all neighbor
    ! pairs only once
end do

do i_cell = 1, n_cells_tot
    call make_cell_list(i_cell, n_cells, n_dim, n_nei_tot, cell_neigh_ls_tot(i_cell,:))
    ! All the neighbors are given, in order to loop over all neighbor
    ! pair
    
    !DEBUG
    !Print*,"inside suboutine get_cell_neigh_list"
    !Print*,"cell_neigh_ls_tot of cell ",i_cell
    !Print*,cell_neigh_ls_tot(i_cell,:)
    !Print*,""
end do


end subroutine get_cell_neigh_list


subroutine make_cell_list(in_cell, n_cells, n_dim, n_nei_tot, cell_list)
#include 'control_simulation.h'
!This routine gives an array with the index of all the neighbouring cells of the
!input cell. Already corrected fo boundary conditions, assuming PBC in all 
!directions
!
implicit none
integer, intent(in) :: in_cell, n_dim, n_cells(n_dim), n_nei_tot
integer, intent(out) :: cell_list(n_nei_tot)
integer :: i_dim, j_dim, j_cell, i_list, x_cell, y_cell, z_cell, i, j, k, l

!!DEBUG
!print*,""
!print*,"Building cell_neigh_ls_tot"
!print*,"in_cell=",in_cell
!print*,"n_cells",n_cells
!print*,"n_nei_tot",n_nei_tot
!!/DEBUG

    z_cell = int( (in_cell - 1) / (n_cells(1)*n_cells(2)) + 1 )
    y_cell = int( ( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - 1 ) / n_cells(1) + 1 )
    x_cell = int( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - n_cells(1) * (y_cell - 1)) 
    i_list = 1
!DEBUG
!print*,"x_cell,y_cell,z_cell",x_cell,y_cell,z_cell    
!DEBUG

select case (n_nei_tot)
    case(27) ! l_cell = r_cut ! 27 neighbors per cell
    do i = 1, 3 ! Go only once over pairs 
       do j = 1, 3 
           do k = 1, 3

                j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) * n_cells(1)

                !DEBUG!
                !print*,"j_cell",j_cell
                !/DEBUG

                !HERE CHECK PBC
                if ( (x_cell .eq. 1) .and. (i.eq.1) ) j_cell = j_cell + n_cells(1)
                if ( (x_cell .eq. n_cells(1)) .and. (i.eq.3) ) j_cell = j_cell - n_cells(1)
                if ( (y_cell .eq. 1) .and. (j.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)
                if ( (y_cell .eq. n_cells(2) ) .and. (j.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2) 
#       if SYMMETRY == 1
                if ( (z_cell .eq. 1) .and. (k.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)*n_cells(3)
                if ( (z_cell .eq. n_cells(3) ) .and. (k.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2)*n_cells(3)   
#       elif  SYMMETRY == 0
                if ( (z_cell + k-2) .lt. 1 ) j_cell = 0 !Don't apply PBC in z, put cell 0 as neighbor
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = 0 !cell 0 is empty, and will be ignored in the 
                                                                  !neighbor cells loop                            
#       endif    

!!DEBUG
!print*,"i_list",i_list
!print*,"i, j, k", i, j, k 
!print*,"cell_list(i_list)", j_cell
!!/DEBUG
                cell_list(i_list) = j_cell
                i_list = i_list + 1
                
            end do
        end do
    end do

    !DEBUG
!print*,"inside subroutine make_cell_list"
!print*,"Cell neigh list of cell", in_cell
!print*,"shape",shape(cell_list)
!print*,"size",size(cell_list)
!print*,"n_nei_tot",n_nei_tot
!print*,cell_list(:)
!print*,""
!/DEBUG

    case(125) ! l_cell = r_cut / 2 ! 62 neighbors per cell
    do i = 0, 4 ! if l_cell = r_cut / 2
        do j = 0, 4 
           do k = 0, 4
               j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) * n_cells(1)

               !HERE CHECK PBC
                !if ( (x_cell+i-2) .lt. 1 ) j_cell = j_cell + n_cells(1)
                if ( (x_cell+i-2) .gt. n_cells(1) ) j_cell = j_cell - n_cells(1)
                if ( (y_cell+j-2) .lt. 1 ) j_cell = j_cell + n_cells(1)*n_cells(2)  
                if ( (y_cell+j-2) .gt. n_cells(2) ) j_cell = j_cell - n_cells(1)*n_cells(2)  
#       if SYMMETRY == 1
                if ( (z_cell + k-2) .lt. 1 ) j_cell = j_cell + n_cells(1)*n_cells(2)*n_cells(3)
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = j_cell - n_cells(1)*n_cells(2)*n_cells(3)   
#       elif  SYMMETRY == 0 
                if ( (z_cell + k-2) .lt. 1 ) j_cell = 0 !Don't apply PBC in z, put cell 0 as neighbor
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = 0   !cell 0 is empty, and will be ignored in the 
                                                                    !neighbor cells loop                            

#       endif                    
                cell_list(i_list) = j_cell
                i_list = i_list + 1
            end do
        end do
    end do

   
    case default
    print*,"Error in routine make_cell_list2, in cell_list"
    stop
    end select
end subroutine
     



subroutine make_cell_list2(in_cell, n_cells, n_dim, n_nei, cell_list)
#include 'control_simulation.h'
!This routine gives an array with the index of all the neighbouring cells of the
!input cell. Already corrected fo boundary conditions, assuming PBC in all 
!directions
!
! OPTIMIZATION DONE FOR n_dim = 3: Loop between neighboring cells in one direction only,
! to count interactions only once. ( do i=2,3; do j=2,3; bla; end do; end do)
! instead of 3**2-1 neighbors, (3**2 -1 )/2 neighbors.
implicit none
integer, intent(in) :: in_cell, n_dim, n_cells(n_dim), n_nei
integer, intent(out) :: cell_list(n_nei)
integer :: i_dim, j_dim, i_cell, j_cell, i_list, x_cell, y_cell, z_cell, i, j, k, l

select case(n_dim)

case(1)
    i_list = 1
    do i = 2, 3
        j_cell = in_cell + i - 2
        if ( j_cell .eq. in_cell ) cycle ! Dont count in_cell as neighbour
    
        !HERE CHECK PBC
        if ( j_cell .eq. 0 ) j_cell = n_cells(1)
        if ( j_cell .eq. ( n_cells(1) + 1 ) ) j_cell = 1

        cell_list(i_list) = j_cell
        i_list = i_list + 1
    end do
    
     
case(2)
    y_cell = (in_cell - 1) / n_cells(1) + 1  ! Get y coordinate of in_cell
    x_cell = in_cell - n_cells(1) * (y_cell - 1)
    i_list = 1
    do i = 2, 3
        do j = 1, 3
            if((i.eq.2).and.(j.eq.3)) cycle
            j_cell = in_cell + i-2 + (j-2) * n_cells(1)
            if ( j_cell .eq. in_cell ) cycle ! don't count same cell

            !HERE CHECK FOR PBC IN ALL DIRECTIONS 
            if ( (x_cell .eq. n_cells(1)) .and. (i.eq.3) ) j_cell = j_cell - n_cells(1)
            if ( (x_cell .eq. 1) .and. (i.eq.1) ) j_cell = j_cell + n_cells(1)
            if ( (y_cell .eq. n_cells(2) ) .and. (j.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2)    
            if ( (y_cell .eq. 1) .and. (j.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)           
            
            cell_list(i_list) = j_cell
            i_list = i_list + 1
        end do
    end do
    
case(3)

    z_cell = int( (in_cell - 1) / (n_cells(1)*n_cells(2)) + 1 )
    y_cell = int( ( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - 1 ) / n_cells(1) + 1 )
    x_cell = int( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - n_cells(1) * (y_cell - 1)) 
    i_list = 1

    select case (n_nei)
    case(13) ! l_cell = r_cut ! 13 neighbors per cell
    do i = 2, 3 ! Go only once over pairs 
       do j = 1, 3 
           if( (i.eq.2) .and. (j.le.1))  cycle !  Go only once over pairs
            do k = 1, 3
               if( (i.eq.2) .and. (j.eq.2) .and. (k.ge.3))  cycle !  Go only once over pairs
                j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) * n_cells(1)
                if ( j_cell .eq. in_cell ) cycle

                !HERE CHECK PBC
                if ( (x_cell .eq. 1) .and. (i.eq.1) ) j_cell = j_cell + n_cells(1)
                if ( (x_cell .eq. n_cells(1)) .and. (i.eq.3) ) j_cell = j_cell - n_cells(1)
                if ( (y_cell .eq. 1) .and. (j.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)
                if ( (y_cell .eq. n_cells(2) ) .and. (j.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2) 
#       if SYMMETRY == 1
                if ( (z_cell .eq. 1) .and. (k.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)*n_cells(3)
                if ( (z_cell .eq. n_cells(3) ) .and. (k.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2)*n_cells(3)   
#       elif  SYMMETRY == 0
                if ( (z_cell + k-2) .lt. 1 ) j_cell = 0 !Don't apply PBC in z, put cell 0 as neighbor
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = 0 !cell 0 is empty, and will be ignored in the 
                                                                  !neighbor cells loop                            
#       endif      

                cell_list(i_list) = j_cell
                i_list = i_list + 1
            end do
        end do
    end do

    case(62) ! l_cell = r_cut / 2 ! 62 neighbors per cell
    do i = 2, 4 ! if l_cell = r_cut / 2
        do j = 0, 4 
            if( (i.eq.2) .and. (j.le.1))  cycle !  Go only once over pairs
            do k = 0, 4
                if( (i.eq.2) .and. (j.eq.2) .and. (k.ge.3))  cycle !  Go only once over pairs
                j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) * n_cells(1)

                if ( j_cell .eq. in_cell ) cycle

                !HERE CHECK PBC
                !if ( (x_cell+i-2) .lt. 1 ) j_cell = j_cell + n_cells(1)
                if ( (x_cell+i-2) .gt. n_cells(1) ) j_cell = j_cell - n_cells(1)
                if ( (y_cell+j-2) .lt. 1 ) j_cell = j_cell + n_cells(1)*n_cells(2)  
                if ( (y_cell+j-2) .gt. n_cells(2) ) j_cell = j_cell - n_cells(1)*n_cells(2)  
#       if SYMMETRY == 1
                if ( (z_cell + k-2) .lt. 1 ) j_cell = j_cell + n_cells(1)*n_cells(2)*n_cells(3)
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = j_cell - n_cells(1)*n_cells(2)*n_cells(3)   
#       elif  SYMMETRY == 0 
                if ( (z_cell + k-2) .lt. 1 ) j_cell = 0 !Don't apply PBC in z, put cell 0 as neighbor
                if ( (z_cell + k-2) .gt. n_cells(3) )  j_cell = 0   !cell 0 is empty, and will be ignored in the 
                                                                    !neighbor cells loop                            

#       endif                    
                cell_list(i_list) = j_cell
                i_list = i_list + 1
            end do
        end do
    end do

   
    case default
    print*,"Error in routine make_cell_list2, in cell_list"
    stop
    end select

case(4)
    i_list = 1
    do i = 1, 3
        do j = 1, 3 
            do k = 1, 3
                do l = 1, 3
                    j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) + (l-2) * n_cells(3)
                    if ( j_cell .eq. in_cell ) cycle
                    !HERE CHECK PBC
                    cell_list(i_list) = j_cell
                    i_list = i_list + 1
                end do
            end do
        end do
    end do

case default
    print*, "ERROR in routine make_cell_list. "
    print*, "This routine works fine for 3 dimensions or less "
    stop

end select 

end subroutine
                     
subroutine get_r_cell(in_cell, n_cells, n_dim, l_cell, r_cell)
!Get the position of the cell in_cell
implicit none
integer, intent(in) :: in_cell, n_dim, n_cells(n_dim)
real(kind=8) , intent(in) :: l_cell(n_dim)
real(kind=8) , intent(out) :: r_cell(n_dim)
integer :: x_cell, y_cell, z_cell

    z_cell = int( (in_cell - 1) / (n_cells(1)*n_cells(2)) + 1 )
    y_cell = int( ( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - 1 ) / n_cells(1) + 1 )
    x_cell = int( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - n_cells(1) * (y_cell - 1)) 

    r_cell(1) = ( float(x_cell) + .5 ) * l_cell(1)
    r_cell(2) = ( float(y_cell) + .5 ) * l_cell(2)
    r_cell(3) = ( float(z_cell) + .5 ) * l_cell(3)

end subroutine

         

!subroutine make_neig_ls(part_in_cell, r_nei, r0, n_part, cell_neigh_ls, n_cells_tot, n_dim, n_cells, l_cell, ff_list)
!!This routine is supposed to do do a neighbor Verlet list, from the cell list. 
!! It has a major disadvantage: the number of neighboring particles for each
!! particle is not known a priori, so a lot of space is wasted
!implicit none
!real(kind=8), intent(in) :: r0(n_dim,n_part), l_cell(n_dim)
!integer, intent(in) :: n_part, n_cells_tot, n_dim, n_cells(n_dim), cell_neigh_ls( n_cells_tot, int((3**n_dim-1)/2)), r_nei(n_part), part_in_cell(n_cells_tot)
!integer, intent(out) :: ff_list(0:256,n_part) ! HERE IS THE PROBLEM: Lot of memory wasted 
!integer :: i_part, j_part, i_cell, j_cell,j
!ff_list(0,:) = 0 ! set number of neighbors to zero for all particles
!!loop over all particles
!do i_part=1, n_part 
!    call get_cell(r0(:,i_part), n_dim, n_cells, l_cell, i_cell) ! get cell of particle i_part
!    !loop over neighboring cells
!    do j=1,3**n_dim-1                                               
!        j_cell = cell_neigh_ls(i_cell,j)                               
!        j_part = part_in_cell(j_cell)                               
!        do while(j_part .ne. 0)     !loop over particles in neighbor cell
!            !if distance between particles is less than cut-off
!            ff_list(0,i_part) = ff_list(0,i_part) + 1
!            ff_list(ff_list(0,i_part),i_part) = j_part         
!            !end if 
!            j_part = r_nei(j_part) !gets next particle in this cell 
!        end do                                                      
!    end do                                                          
!
!! loop over particles in the same cell also (i_cell)           
!
!    j_part = part_in_cell(i_cell)                                   
!    do while(j_part .ne. 0)                                         
!        if(j_part .gt. i_part) then ! count interactions only once 
!            !if distance between particles is less than cut-off
!            ff_list(0,i_part) = ff_list(0,i_part) + 1
!            ff_list(ff_list(0,i_part),i_part) = j_part
!            !end if 
!        end if                                                      
!    
!       j_part = r_nei(j_part) !gets next particle in this cell     
!    end do                                                         
!end do
!
!end subroutine ! make_neig_ls


! neighbor list has the disadvantage of wasting space, which cell list does not
! have. 


!To calculate L-J interactions:
! Variables that should be global to acces: 
! n_cells_tot, part_in_cell, cell_neigh_ls, r_nei, l_nei, n_cells, l_cell
! First loop over cells, because this is parallelizable with OpenMP
! While loops cannot be parallelized
! do i_cell = 1, n_cells_tot ! loop over cells
!     i_part = part_in_cell(i_cell) ! 
!     do while(i_part .ne. 0)  ! loop over particles in this cell
!         do j_dummy = 1, int((3**n_dim-1)/2) ! loop over neighboring cells
!             j_cell = cell_neigh_ls(i_cell,j_dummy)
!             j_part = part_in_cell(j_cell) 
!             do while(j_part .ne. 0)  ! loop over particles in this j_cell  
!                 !Calculate force r0(:,i_part) - r0(:,j_part)
!                 j_part = r_nei(j_part) 
!             end do ! particles in j_cell 
!         end do  ! loop over neighbor cells
!
! Now check interactions in the same cell                            
!         j_part = part_in_cell(i_cell) 
!         do while(j_part .ne. 0)  ! loop over particles in this cell
!             if (i_part .lt. j_part) then !Count interactions just once
!                 !Calculate force r0(:,i_part) - r0(:,j_part)
!             end if    
!             j_part = r_nei(j_part) !gets next particle in this cell     
!         end do  ! loop over neighbor particles in i_cell, j_part
!     end do ! loop over particles i_part in i_cell
! end do !loop over cells                                                    
!











