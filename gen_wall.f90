
  subroutine gen_wall()
#include 'control_simulation.h'
      use commons ;use util ; implicit none
      integer ::info

! ----- Solid wall: Einstein's FCC  
!
!  define wall atoms always to be of type "n_type"
!

          do i_part = n_mon_tot+1,n_part
          a_type(i_part) = n_type
          end do

! ---- Begin generation of wall atoms

!
!---  define equilibrium and start up sites of wall atoms
!
!        info=1
!        if(info.eq.1) then
            i_wall = 1
            i_cell_w_x = n_cell_w_x
            i_cell_w_y = n_cell_w_y
            if((i_cell_w_x.eq.i_cell_w_y).or.                                &
                (i_cell_w_x.eq.2*i_cell_w_y)) then
!----   commensurate surfaces
            do i_cell_w_x = 0,n_cell_w_x-1
                do i_cell_w_y = 0,n_cell_w_y-1
!-----   x coordinates
                    r_wall_equi(1,i_wall) = i_cell_w_x*x_space
                    r_wall_equi(1,i_wall+1) = r_wall_equi(1,i_wall) + x_space/2
                    r_wall_equi(1,i_wall+n_wall/2) = r_wall_equi(1,i_wall)
                    r_wall_equi(1,i_wall+1+n_wall/2) = r_wall_equi(1,i_wall+1)
!----    y coordinates
                    r_wall_equi(2,i_wall) = i_cell_w_y*y_space
                    r_wall_equi(2,i_wall+1) = r_wall_equi(2,i_wall) + y_space/2
                    r_wall_equi(2,i_wall+n_wall/2) = r_wall_equi(2,i_wall)
                    r_wall_equi(2,i_wall+1+n_wall/2) = r_wall_equi(2,i_wall+1)
                    !        call r250(mz,random,n_part,2,kptr)
                    !        r_wall_equi(2,i_wall+n_wall/2) = boundary(2)*random(1)
                    !        r_wall_equi(2,i_wall+1+n_wall/2) = boundary(2)*random(2)
                    i_wall = i_wall + 2
                end do
            end do
        else if(((i_cell_w_x.eq.31).and.(i_cell_w_y.eq.18)).or.          &
            ((i_cell_w_x.eq.62).and.(i_cell_w_y.eq.36))) then
        !----   incommensurate surfaces
        do i_cell_w_x = 0,n_cell_w_x-1
            do i_cell_w_y = 0,n_cell_w_y-1
                !-----   x coordinates
                r_wall_equi(1,i_wall) = i_cell_w_x*x_space
                r_wall_equi(1,i_wall+1) = r_wall_equi(1,i_wall) + x_space/2
                r_wall_equi(2,i_wall+n_wall/2) = r_wall_equi(1,i_wall)
                r_wall_equi(2,i_wall+1+n_wall/2) = r_wall_equi(1,i_wall+1)
                !----    y coordinates
                r_wall_equi(2,i_wall) = i_cell_w_y*y_space
                r_wall_equi(2,i_wall+1) = r_wall_equi(2,i_wall) + y_space/2
                r_wall_equi(1,i_wall+n_wall/2) = r_wall_equi(2,i_wall)
                r_wall_equi(1,i_wall+1+n_wall/2) = r_wall_equi(2,i_wall+1)
                i_wall = i_wall + 2
            end do
        end do
    end if
!   Z coordinates

    do i_wall = 1,n_wall/2
        r_wall_equi(n_dim,i_wall) = z_space_wall
        r_wall_equi(n_dim,i_wall+n_wall/2) = 0.
    end do
      
!!!     else
!!! 
!!!     open(32,file="conf_wall",status="old")
!!! 
!!!     read(32,*) i_wall
!!!     if(i_wall.ne.n_wall) then
!!!         write(*,*) "wrong number of wall atoms in conf_wall"
!!!         stop
!!!     else
!!!         read(32,*)
!!!         do i_wall = 1,n_wall
!!!             i_part = i_wall + n_mon_tot
!!!             read(32,*) a_type(i_part),r_wall_equi(3,i_wall),                 &
!!!                 r_wall_equi(1,i_wall),r_wall_equi(2,i_wall)
!!!         end do
!!!     end if
!!!     close(32)
!!! end if

!
!---  Write equilibrium positions onto actual positions
!

        do i_dim = 1,n_dim
            do i_wall = 1,n_wall
                i_part = n_mon_tot+i_wall
                r0(i_dim,i_part) = r_wall_equi(i_dim,i_wall)
            end do
        end do

  end subroutine gen_wall
