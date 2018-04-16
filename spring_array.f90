subroutine spring_array(mode)
#include "control_simulation.h"
#ifdef SPRING_ARRAY
    ! Routine that introduces an elastic (Hook law) force between neighboring chains
    ! The force is applied to the second bead of the chain
    use commons
    implicit none

    INTEGER,  intent(in) :: mode

    INTEGER, SAVE :: cols,rows, bead_site
    INTEGER :: i,j,l,i_ch,ni,nj,ir0,igraft,ineigh
    REAL (KIND=8) ::  r_neigh(3), r_graft(3), pot_fue(0:3)!,d, dist
    REAL (KIND=8), ALLOCATABLE, SAVE :: f_array(:,:)

    select case (mode)

    case(0)

        cols= NINT(SQRT(0.5*n_chain*boundary(1)/boundary(2)))
        rows= NINT(SQRT(0.5*n_chain*boundary(2)/boundary(1)))

        if ( (rows.eq.1).or.(cols.eq.1) ) then
            print *, " **  WARNING: Either rows or columns are equal to 1. This will break periodic boundary conditions for the spring array "
            print *, "rows=",rows,"  cols=",cols
            print *, " ** Exiting now..."
            STOP
        endif

!        k_spr_x= 800
!        k_spr_y= 0
        ALLOCATE(f_array(3,n_part))
        bead_site= 1 ! This is the number of the monomer in each chain, with grafted =0. It should be ranged
        ! between 1 and n_mon-1. If not, weird array springs will be formed.

        ALLOCATE(Mindex(0:rows+1,0:cols+1))

        Mindex=-1 ! Initializing to -1.

        ! This matrix should contain all indexes of the 1st bead of every chain
        do j=1,cols
            do i=1,rows
                i_ch = 1-i + j*rows 
                Mindex(i,j) = n_mon*(i_ch-1) +1 ! +1 is the grafted bead
            end do
        end do

        ! This is the trick for periodic boundary conditions
        Mindex(:,0) = Mindex(:,cols)
        Mindex(:,cols+1) = Mindex(:,1)
        Mindex(0,:) = Mindex(rows,:)
        Mindex(rows+1,:) = Mindex(1,:)
        Mindex(rows+1,cols+1)=0 
        Mindex(0,0)=0 
        Mindex(rows+1,0)=0 
        Mindex(0,cols+1)=0 
        
        print *, ""
        print *, " * (Spring Array ON) Chains are linked by springs"
        print *, " -- k_spr_x=", k_spr_x 
        print *, " -- k_spr_y=", k_spr_y 
        print *, " -- Index Array:"
        do l=0,rows+1
            print *, Mindex(l,:)
        end do
        print *, " * Lattice parameter", sqrt( 2 * boundary(1) * boundary(2) / n_chain )

    case(1)
        f_array=0.d0
        v_array=0.d0
        do l=0,1 ! Bot/Top loop
            do j=1,cols
                do i=1,rows
                    igraft= Mindex(i,j)+l*n_mon*n_chain/2 
                    ir0= igraft + bead_site

                    do nj=-1,1,2
                        ineigh= Mindex(i,j+nj)+ bead_site +l*n_mon*n_chain/2

                        pot_fue=func_pot_fue(k_spr_x)
                        !!
                        r_neigh = r0(:, ineigh) - r0(:, ir0) 
                        r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 

                        !print *, 'A',' x', r_graft(1), boundary(1)*INT(2.*r_graft(1)*inv_boundary(1))
                        r_graft(1) = r_graft(1) - boundary(1)*INT(2.*r_graft(1)*inv_boundary(1))
                        !print *, 'D',' x', r_graft(1)
                        !print *, 'A',' y', r_graft(2), boundary(2)*INT(2.*r_graft(2)*inv_boundary(2))
                        r_graft(2) = r_graft(2) - boundary(2)*INT(2.*r_graft(2)*inv_boundary(2))
                        !print *, 'D',' y', r_graft(2)
                        
                        !print *, 'A',' x', r_neigh(1), boundary(1)*INT(2.*r_neigh(1)*inv_boundary(1))
                        r_neigh(1) = r_neigh(1) - boundary(1)*INT(2.*r_neigh(1)*inv_boundary(1))
                        !print *, 'D',' x', r_neigh(1)
                        r_neigh(2) = r_neigh(2) - boundary(2)*INT(2.*r_neigh(2)*inv_boundary(2))

                        !print *, NORM2(r_neigh) - NORM2(r_graft)


                        !!

                        !print *, '#', NORM2(r_neigh), igraft
                        !print *, '!',r_neigh

                        f_array(:, ir0) = f_array(:, ir0) + pot_fue(1:3)
                        v_array = v_array + pot_fue(0)
!print *, pot_fue(1:3),'--++'
                    end do !nj

                    do ni=-1,1,2
                        ineigh= Mindex(i+ni,j)+ bead_site +l*n_mon*n_chain/2

                        pot_fue=func_pot_fue(k_spr_y)
                        !!
                        r_neigh = r0(:, ineigh) - r0(:, ir0) 
                        r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
                        r_graft(1) = r_graft(1) - boundary(1)*INT(2.*r_graft(1)*inv_boundary(1))
                        r_graft(2) = r_graft(2) - boundary(2)*INT(2.*r_graft(2)*inv_boundary(2))
                        r_neigh(1) = r_neigh(1) - boundary(1)*INT(2.*r_neigh(1)*inv_boundary(1))
                        r_neigh(2) = r_neigh(2) - boundary(2)*INT(2.*r_neigh(2)*inv_boundary(2))
                        !!
                        !print *, '#', NORM2(r_graft), igraft
                        
                        f_array(:, ir0) = f_array(:, ir0) + pot_fue(1:3)
                        v_array = v_array + pot_fue(0)

!print *, pot_fue(1:3),'--++'
                    end do !ni
!print *, f_array,'---+'
                   ! do nj=-1,1,2
                   !     ineigh= Mindex(i,j+nj)+ bead_site +l*n_mon*n_chain/2
                   !     r_neigh = r0(:, ir0) - r0(:, ineigh)
                   !     ! r_neigh is the vector that joins the 2nd beads of chains

                   !     r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
                   !     ! r_graft is the vector that joins the 1st beads of chains

                   !     d = SQRT(DOT_PRODUCT(r_graft,r_graft))
                   !     ! d es la magnitud de r_graft. Será la distancia de equil.

                   !     dist = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
                   !     ! dist is the distance between the 2nd beads of neighbors

                   !     f_array(:, ir0) = f_array(:, ir0) - k_spr_x*(dist-d)*r_neigh/dist
                   !     v_array = v_array + .5*k_spr_x*(dist-d)**2
                   !     !print *, 'dist=',dist,'||  d=', d
                   !     !print *, k_spr_x*(dist-d)*r_neigh/dist

                   ! end do !nj

                   ! do ni=-1,1,2
                   !     ineigh= Mindex(i,j+ni)+ bead_site +l*n_mon*n_chain/2
                   !     r_neigh = r0(:, ir0) - r0(:, ineigh)
                   !     ! r_neigh is the vector that joins the 2nd beads of chains

                   !     r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
                   !     ! r_graft is the vector that joins the 1st beads of chains

                   !     d = SQRT(DOT_PRODUCT(r_graft,r_graft))
                   !     ! d es la magnitud de r_graft. Será la distancia de equil.

                   !     dist = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
                   !     ! dist is the distance between the 2nd beads of neighbors

                   !     f_array(:, ir0) = f_array(:, ir0) - k_spr_y*(dist-d)*r_neigh/dist
                   !     v_array = v_array + .5*k_spr_y*(dist-d)**2
                   !     !print *, 'dist=',dist,'||  d=', d
                   !     !print *, k_spr_y*(dist-d)*r_neigh/dist
                   ! end do !ni
                end do !i
            end do !j
!print *, f_array
        end do !l
!print *, f_array(:,Mindex(2,1)),'|' ,f_array(:,Mindex(2,2)),'|',f_array(:,Mindex(2,3))
!do i=0,rows+1
!print *, i,')',Mindex(i,:)
!enddo
        force = force + f_array ! f_array is a vector FULL of zeros, may be expensive

        ! NEIGHBOURS:
        ! Each bead has 4 neigh. . According to gen_brush case(5) they're 
        ! generated first in rows and then in cols (row -> fast, col -> slow),
        ! starting from the lower left origin towards up and right.
        ! This means that each (closest) neighboring bead is of the form:
        ! -- Horizontal) ( i_chain +/- n_rows )*n_mon
        ! -- Vertical)   ( i_chain +/- 1 )*n_mon
    end select

CONTAINS

FUNCTION func_pot_fue(k) 

    REAL (kind=8), INTENT (in) :: k
    REAL (kind=8) :: r_neigh(3), r_graft(3), d, dg
    REAL (kind=8), DIMENSION(0:3) :: func_pot_fue

    r_neigh = r0(:, ineigh) - r0(:, ir0) 
    ! Points the nieghbor!!
    ! r_neigh is the vector that joins the 2nd beads of chains

    r_graft = r0(:, igraft) - r0(:, ineigh-bead_site) 
    ! r_graft is the vector that joins the 1st beads of chains
    
    !! Periodic boundary conditions
    !! x)
    r_graft(1) = r_graft(1) - boundary(1)*INT(2.*r_graft(1)*inv_boundary(1))
    !! y)
    r_graft(2) = r_graft(2) - boundary(2)*INT(2.*r_graft(2)*inv_boundary(2))
    !! x)
    r_neigh(1) = r_neigh(1) - boundary(1)*INT(2.*r_neigh(1)*inv_boundary(1))
    !! y)
    r_neigh(2) = r_neigh(2) - boundary(2)*INT(2.*r_neigh(2)*inv_boundary(2))


    dg = SQRT(DOT_PRODUCT(r_graft,r_graft))
    ! dg es la magnitud de r_graft. Será la distancia de equil.

    d = SQRT(DOT_PRODUCT(r_neigh,r_neigh))
    ! d is the distance between the 2nd beads of neighbors

    func_pot_fue(1:3) =  k*(d-dg)*r_neigh/d
    func_pot_fue(0) = .5*k*(d-dg)**2

    !print *, dg, d 

END FUNCTION func_pot_fue

#endif
end subroutine spring_array
