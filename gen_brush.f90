  subroutine gen_brush(mode)
      use commons ;use util 
      use ziggurat, only: rnor,uni
      implicit none
      integer, intent(in) :: mode
      integer :: err_ov, i_ov, j_chain,ini_chain,i_col,i_row, n_col,n_row!, j_part!Add by Kevo to avoid overlap between grafting points
      real :: dis_ov_2, min_dis_2, a_br,a_br_2!sigma(1,1)**2 !Add by Kevo to avoid overlap between grafting points
      real (kind=8) :: fi=0., alpha_ini=.52 !0.78
      real (kind=8) :: xi_ini=0.2, xi_end=2.9, rangle
      !  11º -> 0.2
      ! 166º -> 2.9
!a_br = lattice parameter for brush grafting points
      min_dis_2 = sigma(1,1)**2 ! Sets overlap radius 
!---  Specify position of first monomer in chain
#include "control_simulation.h" 
  select case (mode)

  case(1)  ! ------ Brush in top and bottom walls: CHANNEL

! NOTE: This works with brushes only. Free polymers must be generated in another
! context       

    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        do i_dim = 1,n_dim
            if(i_dim.eq.n_dim) then             ! z coor of polymer head
                if(i_chain.le.n_chain/2) then !top wall
                    r0(i_dim,i_part) = z_space_wall-z_head     !ori2.**(1./6.)
                end if
                if(i_chain.gt.n_chain/2) then !bottom wall
                    r0(i_dim,i_part) = z_head !ori2.**(1./6.)
                end if
            end if

            if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                r0(i_dim,i_part) = random(i_dim)*boundary(i_dim)
            end if
        end do

        !----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
500   continue
            call r250(mz,random,n_part,n_dim,kptr)

            do i_dim = 1,n_dim
                if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                    r0(i_dim,i_part) = r0(i_dim,i_part-1) +    (2*random(i_dim)-1.)*r_chain/sqrt3
                    ! PBC in the plane          
                    if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                    else if(r0(i_dim,i_part).le.0.) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                    end if
                end if !end if x or y 
                if(i_dim.eq.n_dim) then             ! if z coor 

                    if(i_chain.le.n_chain/2) then  ! top wall
#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.5
#           endif                

                    end if  ! end z_coor 
                    if(i_chain.gt.n_chain/2) then !bottom wall

#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#           endif                
     !OBSOLETE                if(solv_flag.eq.0) then
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
     !OBSOLETE                else
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.4
     !OBSOLETE                end if

                    end if
                end if !end z coor
                !      end if !ad_flag

                !!    if(ad_flag.eq.0) then !if adsorbed
                !!             r0(i_dim,i_part) = r0(i_dim,i_part-1) +                       &
                !!        &    (2*random(i_dim)-1.)*r_chain/sqrt3
                !!             if(r0(i_dim,i_part).gt.boundary(i_dim)) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                !!             else if(r0(i_dim,i_part).lt.0) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                !!            end if
                !!        end if !end adsorbed

            end do
            i_dummy = 0
            ! check whether new z coordinates is allright

            if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
            if(i_dummy.eq.1) then
                write(*,*) " New setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
                goto 500
            end if
        end do

    end do  ! ends loop over chains

!debug        call write_conf(1,r0(:,1:n_chain*n_mon),10) ; stop   

case(2)  ! ------ Brush ONLY in bottom wall: DROPLET

    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        do i_dim = 1,n_dim
            if(i_dim.eq.n_dim) then             ! z coor of polymer head
!
! DROPLET: only the bottom wall is populated 
!
                r0(i_dim,i_part) = z_head !ori2.**(1./6.)
            end if
            if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                r0(i_dim,i_part) = random(i_dim)*boundary(i_dim)
            end if
        end do

!----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
            501   continue

        call r250(mz,random,n_part,n_dim,kptr)

        do i_dim = 1,n_dim
            if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
! PBC in the plane          
                    if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                    else if(r0(i_dim,i_part).le.0.) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                    end if
            end if !emd if x or y 
            if(i_dim.eq.n_dim) then             ! if z coor 
#            if SOLVENT == 1 || SOLVENT == 2
                   r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                   r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#            endif
            end if
      !not in drop       end if !end z coor
!          end if
!!obs       if(ad_flag.eq.0) then !if adsorbed
!!obs                r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
!!obs                if(r0(i_dim,i_part).gt.boundary(i_dim)) then
!!obs                 r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
!!obs                else if(r0(i_dim,i_part).lt.0) then
!!obs                 r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
!!obs               end if
!!obs           end if !end adsorbed

        end do
         i_dummy = 0
      ! check whether new z coordinates is allright
        if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
         if(i_dummy.eq.1) then
          write(*,*) "new setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
          goto 501
         end if
        end do
        end do  ! ends loop over chains

  case(3)  ! ------ Brush in top and bottom walls: CHANNEL

! NOTE: This works with brushes only. Free polymers must be generated in another
! context       
!1987
    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        err_ov=1              !Add by Kevo to avoid overlap between grafting points
        i_ov=0
        do while(err_ov.gt.0) !Add by Kevo to avoid overlap between grafting points
            err_ov = 0        !Add by Kevo to avoid overlap between grafting points
            i_ov = i_ov + 1   !Add by Kevo to avoid overlap between grafting points
            do i_dim = 1,n_dim
                if(i_dim.eq.n_dim) then             ! z coor of polymer head
                    if(i_chain.le.n_chain/2) then !top wall
                        r0(i_dim,i_part) = z_space_wall-z_head     !ori2.**(1./6.)
                        ini_chain=1
                    end if
                    if(i_chain.gt.n_chain/2) then !bottom wall
                        r0(i_dim,i_part) = z_head !ori2.**(1./6.)
                        ini_chain = n_chain/2 + 1
                    end if
                end if
                
                if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                    r0(i_dim,i_part) = uni()*boundary(i_dim)
                end if
            end do
            ! CHECK OVERLAP !Add by Kevo to avoid overlap between grafting
            ! points. CHECK ONLY OVERLAP IN SAME WALL!!!!!!
            if(i_ov.ge.5000000) then
                print*, "WARNING: PROGRAM STOPPED!"
                print*, "Error: failed to introduce grafting point number ", i_chain, " after ", i_ov, " attempts"
                stop
                !go to 1987
            end if
            do j_chain = ini_chain,i_chain  !loop over all the chains
                j_part = 1 + (j_chain-1)*n_mon
                dis_ov_2=(r0(1,i_part)-r0(1,j_part))**2+(r0(2,i_part)-r0(2,j_part))**2      
                if(dis_ov_2.le.min_dis_2.and.i_part.ne.j_part) then
                    !print*, "Overlap between grafting points detected"
                    err_ov = err_ov + 1
                    exit
                end if  
            end do !END CHECK OVERLAP !Add by Kevo to avoid overlap between grafting points
        end do !Add by Kevo to avoid overlap between grafting points
        !----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
!500   continue ! OLD AND OBSOLETE
            call r250(mz,random,n_part,n_dim,kptr)

            do i_dim = 1,n_dim
                if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                    r0(i_dim,i_part) = r0(i_dim,i_part-1) +    (2*random(i_dim)-1.)*r_chain/sqrt3
                    ! PBC in the plane          
                    if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                    else if(r0(i_dim,i_part).le.0.) then
                        r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                    end if
                end if !end if x or y 
                if(i_dim.eq.n_dim) then             ! if z coor 

                    if(i_chain.le.n_chain/2) then  ! top wall
#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.4
#           endif                

                    end if  ! end z_coor 
                    if(i_chain.gt.n_chain/2) then !bottom wall

#            if SOLVENT == 1 || SOLVENT == 2 
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.4
#           endif                
     !OBSOLETE                if(solv_flag.eq.0) then
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
     !OBSOLETE                else
     !OBSOLETE                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.4
     !OBSOLETE                end if

                    end if
                end if !end z coor
                !      end if !ad_flag

                !!    if(ad_flag.eq.0) then !if adsorbed
                !!             r0(i_dim,i_part) = r0(i_dim,i_part-1) +                       &
                !!        &    (2*random(i_dim)-1.)*r_chain/sqrt3
                !!             if(r0(i_dim,i_part).gt.boundary(i_dim)) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                !!             else if(r0(i_dim,i_part).lt.0) then
                !!              r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                !!            end if
                !!        end if !end adsorbed

            end do
            i_dummy = 0
            ! check whether new z coordinates is allright

            if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
            if(i_dummy.eq.1) then
                write(*,*) " New setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
            !    goto 500
            end if
        end do

    end do  ! ends loop over chains

 case(4)  ! ------ Brush ONLY in bottom wall: DROPLET
!1984
    do i_chain = 1,n_chain  !loop over all the chains
        i_part = 1 + (i_chain-1)*n_mon
        call r250(mz,random,n_part,n_dim,kptr)
        err_ov=1
        i_ov=0
        do while(err_ov.gt.0)
            err_ov = 0
            i_ov = i_ov + 1
            do i_dim = 1,n_dim
                if(i_dim.eq.n_dim) then             ! z coor of polymer head
    !
    ! DROPLET: only the bottom wall is populated 
    !
                    r0(i_dim,i_part) = z_head !ori2.**(1./6.)
                end if
                if(i_dim.ne.n_dim) then                 ! x y coordinates of heads: at random in the plane
                    r0(i_dim,i_part) = uni()*boundary(i_dim)
                end if
            end do
            ! CHECK OVERLAP !Add by Kevo to avoid overlap between grafting
            ! points. CHECK ONLY OVERLAP IN SAME WALL!!!!!!
            if(i_ov.ge.5000000) then
                print*, "WARNING: PROGRAM STOPPED!"
                print*, "Error: failed to introduce grafting point number ", i_chain, " after ", i_ov, " attempts"
                stop
                !go to 1984
            end if
            do j_chain = ini_chain,i_chain  !loop over all the chains
                j_part = 1 + (j_chain-1)*n_mon
                dis_ov_2=(r0(1,i_part)-r0(1,j_part))**2+(r0(2,i_part)-r0(2,j_part))**2      
                if(dis_ov_2.le.min_dis_2.and.i_part.ne.j_part) then
                    !print*, "Overlap between grafting points detected"
                    err_ov = err_ov + 1
                    exit
                end if  
            end do !END CHECK OVERLAP !Add by Kevo to avoid overlap between grafting points
        end do !Add by Kevo to avoid overlap between grafting points


!----  generate random walk for next monomers on chain

        do i_mon = 1,n_mon-1
            i_part = i_part + 1
            504   continue

            call r250(mz,random,n_part,n_dim,kptr)

            do i_dim = 1,n_dim
                if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                    r0(i_dim,i_part) = r0(i_dim,i_part-1) +(2*uni()-1.)*r_chain/sqrt3
                         !(2*random(i_dim)-1.)*r_chain/sqrt3
                        ! PBC in the plane          
                        if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                        else if(r0(i_dim,i_part).le.0.) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                        end if
                end if !emd if x or y 
                if(i_dim.eq.n_dim) then             ! if z coor 
#                if SOLVENT == 1 || SOLVENT == 2
                       r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#                elif SOLVENT == 0 || SOLVENT == 3
                       r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#                endif
                end if
            end do
             i_dummy = 0
      !     check whether new z coordinates is allright
            if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
             if(i_dummy.eq.1) then
              write(*,*) "new setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
              goto 504
             end if
            end do
        end do  ! ends loop over chains

case(5)  ! ------ Ordered Brush in bottom and top wall

!print*, "entro" !DEBUG

    a_br = sqrt( 2 * boundary(1) * boundary(2) / n_chain ) !lattice parameter
    a_br_2 = a_br / 2

!Check parameter compatibility for ordered brush
    n_col = nint(boundary(1)/a_br)
    if(mod(boundary(1),a_br).ne.0.) then 
        print*, "boundary(1), should be a multiple of the lattice parameter a_br"
        print*, "boundary(1) = ",boundary(1),"a_br = ",a_br
        print*,"boundary(1)/a_br = ", boundary(1)/a_br
    end if

    n_row = nint(boundary(2)/a_br)
    if(mod(boundary(2),a_br).ne.0.) then
       print*, "boundary(2), should be a multiple of the lattice parameter a_br"
       print*, "boundary(2) = ",boundary(2),"a_br = ",a_br
       print*,"boundary(2)/a_br = ", boundary(2)/a_br
    end if
    
    if(n_chain.ne.2*n_col*n_row) then
        print*,"ERROR: n_chain should be equal to 2*n_col*n_row"
        print*,"n_chain: ",n_chain
        print*,"2*n_col*n_row= ","2*",n_col,"*",n_row," = ", 2*n_col*n_row 
        stop
    end if

!Set head positions
do j_chain = 0, 1 !0 is top wall. 1 is for bottom wall
    i_part =  j_chain * n_mon * n_chain / 2 ! starting particle to set it's location
    do i_col = 0, n_col-1 !loop for columns
        do i_row = 0, n_row-1 !loop for rows
            i_part = i_part + 1
            r0(1,i_part) = a_br_2 + a_br*i_col !set head ubication of top grafting point
            r0(2,i_part) = a_br_2 + a_br*i_row 
            if(i_part.le.n_mon*n_chain/2) then 
                r0(3,i_part) = z_space_wall-z_head
            else if(i_part.gt.n_mon*n_chain/2) then
                r0(3,i_part) = z_head
            end if
!            print*,i_part,r0(:,i_part) ! DEBUG
            !Set the remaining monomers in chain
            do i_mon=2,n_mon
                i_part = i_part + 1
                do i_dim=1,3
                    if(i_dim.lt.n_dim) then    !x and y                
                        r0(i_dim,i_part) = r0(i_dim,i_part-1) +(2*uni()-1.)*r_chain/sqrt3 !bottom
                        ! PBC in the plane          
                        if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                        else if(r0(i_dim,i_part).le.0.) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                        end if
                        !print*,"i_dim",i_dim,"n_dim",n_dim !DEBUG
                    else if(i_dim.eq.n_dim) then !z
                         !print*,"i_dim = ",n_dim ! DEBUG
                         if(i_part.le.n_mon*n_chain/2) then  ! top wall
                         !print*,"TOP WALL" !DEBUG
#            if SOLVENT == 1 || SOLVENT == 2 
                             r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.5
#           endif       
                        else if (i_part.gt.n_mon*n_chain/2) then !bottom wall    
                        !print*, "BOTTOM WALL" !DEBUG
#            if SOLVENT == 1 || SOLVENT == 2 
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#           endif                
                        end if  ! end top or bottom 
                    end if ! end x,y or z
                end do ! end i_dim
!                print*,i_part,r0(:,i_part) !DEBUG
            end do  ! end rest of the chain          
        end do !end i_row
    end do ! end i_col
end do ! end j_chain bottom or top wall


!print*,"salio" !DEBUG

case(6)  ! ------ Ordered Brush ONLY in bottom wall: DROPLET
    a_br = sqrt( boundary(1)*boundary(2)/n_chain ) !lattice parameter
    a_br_2 = a_br / 2

!Check parameter compatibility for ordered brush
    if(mod(boundary(1),a_br).eq.0.) then 
        n_col = int(boundary(1)/a_br)
    else
        print*, "ERROR: number of chains incorrect. Stopping programm"
        print*, "boundary(1), should be a multiple of the lattice parameter a_br"
        print*, "boundary(1) = ",boundary(1),"a_br = ",a_br
        print*,"boundary(1)/a_br = ", boundary(1)/a_br
        stop 
    end if

    if(mod(boundary(2),a_br).eq.0.) then
        n_row = int(boundary(2)/a_br)
    else
       print*, "ERROR: number of chains incorrect. Stopping programm"
       print*, "boundary(2), should be a multiple of the lattice parameter a_br"
       print*, "boundary(2) = ",boundary(2),"a_br = ",a_br
       print*,"boundary(2)/a_br = ", boundary(2)/a_br
       stop 
    end if
   
    q_part=1! grafted particle 
    do i_col = 0, n_col-1 !loop for columns
        do i_row = 0, n_row-1 !loop for rows
            i_part = q_part  ! i_part reference to build the rest of the chain
            r0(1,i_part) = a_br_2 + a_br*i_col !set head ubication
            r0(2,i_part) = a_br_2 + a_br*i_row
            
            q_part = q_part + n_mon !next head 
            r0(3,i_part) = z_head

                do i_mon = 1,n_mon-1
                    i_part = i_part + 1
                    509 continue 
                    call r250(mz,random,n_part,n_dim,kptr)

                        do i_dim = 1,n_dim
                            if(i_dim.lt.n_dim) then ! for x and y coordinates, random
                                r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
! PBC in the plane          
                                if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                                    r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                                else if(r0(i_dim,i_part).le.0.) then
                                    r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                                end if
                            end if !emd if x or y 
                            if(i_dim.eq.n_dim) then             ! if z coor 
#            if SOLVENT == 1 || SOLVENT == 2
                               r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
#            elif SOLVENT == 0 || SOLVENT == 3
                               r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
#            endif
                             end if
                        end do
                        i_dummy = 0
                     ! check whether new z coordinates is allright
                       if((r0(n_dim,i_part).le.0.).or.(r0(n_dim,i_part).ge.boundary(n_dim)))  i_dummy = 1
                        if(i_dummy.eq.1) then
                         write(*,*) "new setup attempt for monomer (z-pos)",  i_part,r0(n_dim,i_part)
                         goto 509
                        end if
                       end do
        end do
    end do    

case(7)  ! ------ Ordered Brush in bottom and top wall, starting with initial orientation
! defined by 'alpha_ini' (angle between two consecutive beads and the top/bot wall) and
! 'fi' (angle between the projection of two consecutive beads over top/bot wall and the x axis)

print*, "entro" !DEBUG

    a_br = sqrt( 2 * boundary(1) * boundary(2) / n_chain ) !lattice parameter
    a_br_2 = a_br / 2

    print*, "  * Lattice parameter a_br = ", a_br

!Check parameter compatibility for ordered brush
    n_col = nint(boundary(1)/a_br)
    if(mod(boundary(1),a_br).ne.0.) then 
        print*, "boundary(1), should be a multiple of the lattice parameter a_br"
        print*, "boundary(1) = ",boundary(1),"a_br = ",a_br
        print*,"boundary(1)/a_br = ", boundary(1)/a_br
    end if

    n_row = nint(boundary(2)/a_br)
    if(mod(boundary(2),a_br).ne.0.) then
       print*, "boundary(2), should be a multiple of the lattice parameter a_br"
       print*, "boundary(2) = ",boundary(2),"a_br = ",a_br
       print*,"boundary(2)/a_br = ", boundary(2)/a_br
    end if
    
    if(n_chain.ne.2*n_col*n_row) then
        print*,"ERROR: n_chain should be equal to 2*n_col*n_row"
        print*,"n_chain: ",n_chain
        print*,"2*n_col*n_row= ","2*",n_col,"*",n_row," = ", 2*n_col*n_row 
        stop
    end if

!Set head positions
do j_chain = 0, 1 !0 is top wall. 1 is for bottom wall
    i_part =  j_chain * n_mon * n_chain / 2 ! starting particle to set it's location
    do i_col = 0, n_col-1 !loop for columns
        do i_row = 0, n_row-1 !loop for rows
            i_part = i_part + 1
            r0(1,i_part) = a_br_2 + a_br*i_col !set head ubication of top grafting point
            r0(2,i_part) = a_br_2 + a_br*i_row 
            if(i_part.le.n_mon*n_chain/2) then 
                r0(3,i_part) = z_space_wall-z_head
            else if(i_part.gt.n_mon*n_chain/2) then
                r0(3,i_part) = z_head
            end if
!            print*,i_part,r0(:,i_part) ! DEBUG
            !Set the remaining monomers in chain
#           ifdef RANGLE
            rangle= (xi_end-xi_ini)*uni() + xi_ini
#           else
            rangle=alpha_ini
#           endif
            do i_mon=2,n_mon
                i_part = i_part + 1
                do i_dim=1,3
                    if(i_dim.lt.n_dim) then    !x and y                
                        if(i_dim.eq.1) then    ! x
                                !r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*COS(fi)*COS(alpha_ini)  !bottom
                                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*COS(fi)*COS(rangle)  !bottom
                        else                   ! y 
                                !r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(fi)*COS(alpha_ini) !bottom
                                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(fi)*COS(rangle) !bottom
                        end if
                        ! PBC in the plane          
                        if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                        else if(r0(i_dim,i_part).le.0.) then
                            r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                        end if
                        !print*,"i_dim",i_dim,"n_dim",n_dim !DEBUG
                    else if(i_dim.eq.n_dim) then !z
                         !print*,"i_dim = ",n_dim ! DEBUG
                         if(i_part.le.n_mon*n_chain/2) then  ! top wall
                         !print*,"TOP WALL" !DEBUG
#            if SOLVENT == 1 || SOLVENT == 2 
                             !r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.96*SIN(alpha_ini)
                             r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.96*SIN(rangle)
#            elif SOLVENT == 0 || SOLVENT == 3
                             !r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.96*SIN(alpha_ini)
                             r0(i_dim,i_part) = r0(i_dim,i_part-1) - 0.96*SIN(rangle)
#           endif       
                         else if (i_part.gt.n_mon*n_chain/2) then !bottom wall    
                        !print*, "BOTTOM WALL" !DEBUG
#            if SOLVENT == 1 || SOLVENT == 2 
                            !r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(alpha_ini)
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(rangle)
#            elif SOLVENT == 0 || SOLVENT == 3
                            !r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(alpha_ini)
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(rangle)
#           endif                
                         end if  ! end top or bottom 
                    end if ! end x,y or z
                end do ! end i_dim
!                print*,i_part,r0(:,i_part) !DEBUG
            end do  ! end rest of the chain          
        end do !end i_row
    end do ! end i_col
end do ! end j_chain bottom or top wall

case(8)  ! ------ Ordered Brush ONLY in bottom wall: DROPLET
    a_br = sqrt( boundary(1)*boundary(2)/n_chain ) !lattice parameter
    a_br_2 = a_br / 2

!Check parameter compatibility for ordered brush
    if(mod(boundary(1),a_br).eq.0.) then 
        n_col = int(boundary(1)/a_br)
    else
        print*, "ERROR: number of chains incorrect. Stopping programm"
        print*, "boundary(1), should be a multiple of the lattice parameter a_br"
        print*, "boundary(1) = ",boundary(1),"a_br = ",a_br
        print*,"boundary(1)/a_br = ", boundary(1)/a_br
        stop 
    end if

    if(mod(boundary(2),a_br).eq.0.) then
        n_row = int(boundary(2)/a_br)
    else
       print*, "ERROR: number of chains incorrect. Stopping programm"
       print*, "boundary(2), should be a multiple of the lattice parameter a_br"
       print*, "boundary(2) = ",boundary(2),"a_br = ",a_br
       print*,"boundary(2)/a_br = ", boundary(2)/a_br
       stop 
    end if
   
    q_part=1! grafted particle 
    do i_col = 0, n_col-1 !loop for columns
        do i_row = 0, n_row-1 !loop for rows
            i_part = q_part  ! i_part reference to build the rest of the chain
            r0(1,i_part) = a_br_2 + a_br*i_col !set head ubication
            r0(2,i_part) = a_br_2 + a_br*i_row
            
            q_part = q_part + n_mon !next head 
            r0(3,i_part) = z_head

                do i_mon = 1,n_mon-1
                    i_part = i_part + 1
                        do i_dim = 1,n_dim
                            if(i_dim.eq.1) then ! for x and y coordinates, random
                                !r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
                                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*COS(fi)*COS(alpha_ini)  !bottomi
                            else 
                                r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(fi)*COS(alpha_ini)  !bottomi
! PBC in the plane          
                                if(r0(i_dim,i_part).ge.boundary(i_dim)) then
                                    r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
                                else if(r0(i_dim,i_part).le.0.) then
                                    r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
                                end if
                            end if !emd if x or y 
                            if(i_dim.eq.n_dim) then             ! if z coor 
#            if SOLVENT == 1 || SOLVENT == 2
                             !  r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.2
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(alpha_ini)
#            elif SOLVENT == 0 || SOLVENT == 3
                            !   r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.5
                            r0(i_dim,i_part) = r0(i_dim,i_part-1) + 0.96*SIN(alpha_ini)
#            endif
                             end if
                        end do
                       end do
        end do
    end do    
!print*,"salio" !DEBUG
 end select
end subroutine gen_brush
