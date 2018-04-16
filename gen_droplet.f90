subroutine gen_droplet(mode)
#include 'control_simulation.h'
! This routine is more or less the piece of code
! used in conf_default to generate the first configuration 
! of brushes
! Droplet or melt first conf. generation

        use commons ; use util ; use ziggurat; implicit none
        logical, parameter :: debug = .false.
        real (kind=8) :: r_neigh(3)
        integer, intent(in) :: mode
!the droplet monomers have positions in the interval
! [n_mon*n_chain,nmon_tot]
! i_part value at the begginning of the droplet "segment" of r0

select case (mode)

    case(1) ! generates the droplet or melt drop chains. Particle 3

!        part_init_d = n_mon*n_chain ! the next i_part will be a droplet monomer 
!
!---  Specify position of first monomer in each chain of the droplet

    do i_chain = 1,n_chain_d
       i_part = part_init_d + 1 + (i_chain-1)*n_mon_d
       call r250(mz,random,n_part,n_dim,kptr)
!       do i_dim = 1,n_dim
          r0(:,i_part) = random(:)*boundary(:)
! we leave a skin in z direction          
        if(r0(3,i_part)<= z_skin) r0(3,i_part) = r0(3,i_part) +  z_skin 
        if(r0(3,i_part)>= z_space_wall-z_skin) r0(3,i_part) = r0(3,i_part) -  z_skin 
          !debug
          !print '(3f14.2)',r0(:,i_part)
          !if((mod(i_part-part_init_d,10).eq.0).and.(i_part-part_init_d.gt.10)) print *,' '
!       end do
!----  generate random walk for next monomers on chain
       do i_mon = 1,n_mon_d-1
        i_part = i_part + 1
        i_dummy=0
22     continue
       call r250(mz,random,n_part,n_dim,kptr) !if I understand, random has dimension=3 
!   Adds next monomer with random direction and at a distance of 0.95*R0 of the previous one  
!   but without going inside the skin in z direction 
         r_neigh(:) = (2*random(:)-1.)/sqrt3
         r_neigh(:) = 0.95*r_chain*r_neigh(:)/sqrt(dot_product(r_neigh,r_neigh))
         vec_dummy(:) =  r0(:,i_part-1) + r_neigh(:)
         ! define a counter of trials, to give up if the algor kept trapped 
    i_dummy = i_dummy + 1
 if(i_dummy>1000) then
     print *, "A succesfull conf could not be generated"
     print*,"Stopping here "
     stop
 end if
! new trial is the skin in z is occupied         

         if(vec_dummy(3) <= z_skin .or. vec_dummy(3) > boundary(3)-z_skin) goto 22 
        r0(:,i_part) = vec_dummy(:) 

!   Put the particles inside the boundaries
            do i_dim=1,n_dim
            if(r0(i_dim,i_part).ge.boundary(i_dim)) then
               r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
            else if(r0(i_dim,i_part).le.0.) then
               r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
            end if
            end do

       end do
      end do                !particles loop

      print*,"  *  First conf for melt generated ! "
!debug
!debug if (debug) then
!         call write_conf(1,r0(:,n_mon*n_chain+1:n_mon_tot),n_mon_d) ; stop
!         call write_conf(1,r0(:,:),n_mon_d) ; stop
!debug end if 

   case (2)  ! Generates a configuration of particle 4 type (chains or monomers ) 
! 

! WARN: this is dirty ! redefinition of part_init_d for the different cases.
!obsolete       part_init_d = part_init_e
!
!---  Specify position of first monomer in each chain of the droplet

 do i_chain = 1,n_chain_e
       !ORI: i_part = part_init_d + 1 + (i_chain-1)*n_mon_e
       i_part = part_init_e + 1 + (i_chain-1)*n_mon_e
       call r250(mz,random,n_part,n_dim,kptr)

          r0(:,i_part) = random(:)*boundary(:)

! we leave a skin in z direction          

        if(r0(3,i_part)<= z_skin) r0(3,i_part) = r0(3,i_part) +  z_skin 
        if(r0(3,i_part)>= z_space_wall-z_skin) r0(3,i_part) = r0(3,i_part) -  z_skin 

!----  generate random walk for next monomers on chain

   do i_mon = 1,n_mon_e-1

        i_part = i_part + 1
        i_dummy=0

23     continue

       call r250(mz,random,n_part,n_dim,kptr) !if I understand, random has dimension=3 
!   Adds next monomer with random direction and at a distance of 0.95*R0 of the previous one  
!   but without going inside the skin in z direction 

         r_neigh(:) = (2*random(:)-1.)/sqrt3
         r_neigh(:) = 0.95*r_chain*r_neigh(:)/sqrt(dot_product(r_neigh,r_neigh))
         vec_dummy(:) =  r0(:,i_part-1) + r_neigh(:)
         ! define a counter of trials, to give up if the algor kept trapped 
    i_dummy = i_dummy + 1
 if(i_dummy>1000) then
     print *, "Particle 4: A succesfull conf could not be generated"
     print*,"Stopping here "
     stop
 end if
! new trial is the skin in z is occupied         

         if(vec_dummy(3) <= z_skin .or. vec_dummy(3) > boundary(3)-z_skin) goto 23 
        r0(:,i_part) = vec_dummy(:) 

!   Put the particles inside the boundaries
          do i_dim=1,n_dim
            if(r0(i_dim,i_part).ge.boundary(i_dim)) then
               r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
            else if(r0(i_dim,i_part).le.0.) then
               r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
            end if
          end do

       end do
   end do                !particles loop

      print*,"   * First conf for particle-4 molecules  generated ! "


   case(3)  ! Generates a configuration of particle 2 type (chains or monomers ) for system 2
! 

!obsolete       part_init_d = 0
!
!---  Specify position of first monomer in each chain of the droplet

 do i_chain = 1,n_chain
!ORI: watch the logic !!!!       i_part = part_init_d + 1 + (i_chain-1)*n_mon
       i_part =  1 + (i_chain-1)*n_mon
       call r250(mz,random,n_part,n_dim,kptr)

          r0(:,i_part) = random(:)*boundary(:)

! we leave a skin in z direction          

        if(r0(3,i_part)<= z_skin) r0(3,i_part) = r0(3,i_part) +  z_skin 
        if(r0(3,i_part)>= z_space_wall-z_skin) r0(3,i_part) = r0(3,i_part) -  z_skin 

!----  generate random walk for next monomers on chain

       do i_mon = 1,n_mon-1
        i_part = i_part + 1
        i_dummy=0
24     continue

       call r250(mz,random,n_part,n_dim,kptr) !if I understand, random has dimension=3 
!   Adds next monomer with random direction and at a distance of 0.95*R0 of the previous one  
!   but without going inside the skin in z direction 

         r_neigh(:) = (2*random(:)-1.)/sqrt3
         r_neigh(:) = 0.95*r_chain*r_neigh(:)/sqrt(dot_product(r_neigh,r_neigh))
         vec_dummy(:) =  r0(:,i_part-1) + r_neigh(:)
         ! define a counter of trials, to give up if the algor kept trapped 
    i_dummy = i_dummy + 1
 if(i_dummy>1000) then
     print *, "Particle 2: A succesfull conf could not be generated"
     print*,"Stopping here "
     stop
 end if
! new trial is the skin in z is occupied         

         if(vec_dummy(3) <= z_skin .or. vec_dummy(3) > boundary(3)-z_skin) goto 24
        r0(:,i_part) = vec_dummy(:) 

!   Put the particles inside the boundaries
            do i_dim=1,n_dim
            if(r0(i_dim,i_part).ge.boundary(i_dim)) then
               r0(i_dim,i_part) = r0(i_dim,i_part) - boundary(i_dim)
            else if(r0(i_dim,i_part).le.0.) then
               r0(i_dim,i_part) = r0(i_dim,i_part) + boundary(i_dim)
            end if
            end do

       end do
      end do                !particles loop

      print*,"   * First conf for particle-2 molecules  generated ! "

!debug
        if (debug) then
                call write_conf(1,r0(:,1:part_init_d),n_mon)
                print '(/a/)' , "  *  Writing configuration for particle 2 and quiting  "
                stop
        end if 

    case(4)  ! Droplet generation. Taken from mfa_drop_force
#if SYSTEM == 1

! We put the droplet in the middle of simulation box.

x_droplet_shift = boundary(1)/2. - boundary_d(1)/2.

! check congruence: droplet can not be bigger that the brush layer

if(x_droplet_shift < 0.) then 
print '(/a/)', " Droplet x can not be bigger than polymer brush layer"
stop
end if


! Check that boundary for brush layer and droplet are the same in y direction

    if (boundary(2) /= boundary_d(2)) then
    print '(/a/)' ," The intention is applying PBC in y : dimensions of droplet and layer brush should be equal"
    stop
    end if
! -----------------------------------------
! NOTE: here random has dimension n_part but only the first 3 components are filled with
!       random numbers.
!---------------------------------------------------------------
! NOTE: It uses var. z_min_d for  placing the heads of droplet polymers
!---------------------------------------------------------------


!---  Specify position of first monomer in each chain of the droplet
  do i_chain = 1,n_chain_d
  i_part = part_init_d + 1 + (i_chain-1)*n_mon_d

       call r250(mz,random,n_part,n_dim,kptr)

!       r0(1:2,i_part) = random(1:2)*boundary_d(1:2) 
       r0(1,i_part) = uni()*boundary_d(1)
       r0(2,i_part) = uni()*boundary_d(2)     
       r0(3,i_part) = z_min_d + uni()*( boundary_d(3) - z_min_d )
!DEBUG
 
if((r0(1,i_part).gt.boundary_d(1)).or.(r0(1,i_part).lt.0.)) then
    print*,"lio en x",i_part 
end if

if((r0(2,i_part).gt.boundary_d(2)).or.(r0(2,i_part).lt.0.)) then
    print*,"lio en y",i_part
end if

if((r0(3,i_part).gt.boundary_d(3)).or.(r0(3,i_part).lt.z_min_d)) then
    print*,"lio en z",i_part
end if

!/DEBUG

!
!  leave a skin in z direction          
!
if(r0(3,i_part)<= z_skin) r0(3,i_part) = r0(3,i_part) +  z_skin 
if(r0(3,i_part)>= z_space_wall-z_skin) r0(3,i_part) = r0(3,i_part) -  z_skin 
!
!----  generate random walk for next monomers on chain
!
       do i_mon = 1,n_mon_d-1
        i_part = i_part + 1
        i_dummy=0
25     continue
       call r250(mz,random,n_part,n_dim,kptr) 
       
!   Adds next monomer with random direction and at a distance of 0.95*R0 of the previous one  
!   but without going inside the skin in z direction 

        ! r_neigh(:) = (2*random(1:3)-1.)/sqrt3
        r_neigh(:) = (2*(/uni(),uni(),uni()/)-1.)/sqrt3 
        r_neigh(:) = 0.95*r_chain*r_neigh(:)/sqrt(dot_product(r_neigh,r_neigh))
         vec_dummy(:) =  r0(:,i_part-1) + r_neigh(:)
         ! define a counter of trials, to give up if the algor kept trapped 
         i_dummy = i_dummy + 1

          if(i_dummy>1000) then
              print *, "A succesful conf could not be generated"
              print*,"Stopping here "
!DEB Debug brush configuration

          call write_conf(1,r0(:,:),10)
          stop
          end if
! new trial is the skin in z is occupied         

          if(vec_dummy(3) <= z_skin .or. vec_dummy(3) > boundary(3)-z_skin) goto 25 
          
          r0(:,i_part) = vec_dummy(:) ! if everything is ok, take this position 
        
! DROPLET: only y 
!   Put the particles inside the boundaries
! apply PBC 
            if(r0(2,i_part).ge.boundary_d(2)) then
               r0(2,i_part) = r0(2,i_part) - boundary_d(2)
            else if(r0(2,i_part).le.0.) then
               r0(2,i_part) = r0(2,i_part) + boundary_d(2)
            end if



       end do           ! particles in a chain loop 
  end do                ! chains loop

!
! Shift x of the droplet to the middle of the sample 
!

do i_part=part_init_d +1,n_mon_tot
    r0(1,i_part) = r0(1,i_part) + x_droplet_shift
end do

!DEBUG
!print*,r0
!print*,"Fin posiciones iniciales"
#endif /*system=1; droplet */
#ifdef STARS
case(5)
      do i_part = n_chain*n_mon+n_mon_d*n_chain_d+n_chain_e*n_mon_e+1,n_mon_tot
       if(mod(i_part-(n_chain*n_mon+n_mon_d*n_chain_d+n_chain_e*n_mon_e+1),(n_mon_arm*n_arms+1)).eq.0) then
       a_type(i_part) = 5
       else
       a_type(i_part) = 6
       endif
      end do 


! Generate stars
!




i_part=n_mon*n_chain+n_mon_d*n_chain_d+n_mon_e*n_chain_e


do i_star=1,n_stars

i_part=i_part+1


call r250(mz,random,n_part,n_dim,kptr)

do i_dim = 1,n_dim-1
r0(i_dim,i_part)=random(i_dim)*boundary(i_dim)
r0_star(i_star,i_dim)=r0(i_dim,i_part)
end do


r0(n_dim,i_part)=random(n_dim)*(boundary(n_dim)-2*z_skin)+z_skin
r0_star(i_star,n_dim)=r0(n_dim,i_part)



do i_arm = 1,n_arms
do i_mon_arm = 1, n_mon_arm

i_part=i_part+1


i_dummy=0

600 continue

if (i_dummy.eq.1000) then 
write(*,*) "configuration could not be created"
stop
end if


call r250(mz,random,n_part,n_dim,kptr)


if(i_mon_arm.gt.1) then
do i_dim = 1,n_dim
r0(i_dim,i_part) = r0(i_dim,i_part-1) + (2*random(i_dim)-1.)*r_chain/sqrt3
end do


if(r0(n_dim,i_part).lt.0.5.or.r0(n_dim,i_part).gt.(boundary(n_dim)-0.5)) then
write(*,*) "new setup attempt for star monomer",i_part
i_dummy=i_dummy+1
goto 600
end if



else
do i_dim = 1,n_dim
r0(i_dim,i_part) = r0_star(i_star,i_dim) + (2*random(i_dim)-1.)*r_chain/sqrt3
end do
endif


end do
end do
end do



#endif
end select 


end subroutine gen_droplet
