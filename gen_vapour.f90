subroutine gen_vapour()
! This routine is more or less the piece of code
! used in conf_default to generate the first configuration 
! of brushes
! Low density (or not)  vapour phase generation

       use commons ; use util 
       use functions ! write_conf 
       implicit none
       logical, parameter :: debug = .false.
       integer :: i_part_vap
       real (kind=8) :: r_neigh(3), boundary_vap(3)

! The vapour monomers have positions in the interval
! [n_mon*n_chain+n_chain_f*n_mon_d+1,nmon_tot]

! i_part value at the beginning of the droplet "segment" of r0

          i_part_vap = n_mon*n_chain + n_mon_d*n_chain_d ! the next i_part will be a droplet monomer from here on.
!

! We put the droplet in the middle of the sample of the brushes.
!--------------------------------------
! NOTE: the 'box' for the vapor is taken to be (look for ![***])
!        the space to the left of the 'droplet' and with z_skin to the box.
!        The same is true when the vapor is translated to the right part of the sample.
!--------------------------------------
boundary_vap(1)  = boundary(1)/2. - boundary_d(1)/2. -z_skin ![***]
boundary_vap(2)  = boundary(2)
boundary_vap(3)  = boundary(3) - z_skin


! NOTE: here random has dimension n_part but only the first 3 components are filled with
!       random numbers.
! NOTE: The algorhytm calculates first the first atom of the chain. Then replicates it to the other
!       box of the vapour face. Then the following atoms of the chains are calculated nd replicated. Only
!       for those ones, the applicatin of PBC are necessary. Atom 1 fullfils PBC as defined.


!---  Specify position of first monomer in each chain of the droplet.
! This is done only for the 'left vapour box' and the copied to the right box.

do i_chain = 1,n_chain_vap/2
  
  i_part = i_part_vap + 1 + (i_chain-1)*n_mon_d

       call r250(mz,random,n_part,n_dim,kptr)

       r0(:,i_part) = random(1:3)*boundary_vap(1:3)

!
!  leave a skin in z direction          
!
       if(r0(3,i_part)<= z_skin) r0(3,i_part) = r0(3,i_part) +  z_skin 
       if(r0(3,i_part)>= z_space_wall-z_skin) r0(3,i_part) = r0(3,i_part) -  z_skin 

! First particle of the chain in the right box of vapour

        r0(i_part + n_chain_vap*n_mon_d/2,1) = r0(1,i_part) +  boundary(1)/2. + boundary_d(1)/2. + z_skin ![***]
        r0(i_part + n_chain_vap*n_mon_d/2,2:3) = r0(2:3,i_part)  

!
!----  generate random walk for next monomers on chain
!
    do i_mon = 1,n_mon_d-1
         i_part = i_part + 1
         i_dummy=0
22       continue

         call r250(mz,random,n_part,n_dim,kptr) 
       
!   Adds next monomer with random direction and at a distance of 0.95*R0 of the previous one  
!   but without going inside the skin in z direction 

         r_neigh(:) = (2.*random(1:3)-1.)/sqrt3
         r_neigh(:) = 0.95*r_chain*r_neigh(:)/sqrt(dot_product(r_neigh,r_neigh))
         vec_dummy(:) =  r0(:,i_part-1) + r_neigh(:)
         ! define a counter of trials, to give up if the algor kept trapped 
         i_dummy = i_dummy + 1

         if(i_dummy>1000) then
             print *, "A succesful vapour conf. could not be generated"
             print*,"Stopping here "
             call write_conf(1,r0(:,i_part_vap+1:n_mon_tot),10)
             stop
         end if
!-----  new trial if the skin in z is occupied         

         if(vec_dummy(3) <= z_skin .or. vec_dummy(3) > boundary_vap(3)) goto 22 
          
         r0(:,i_part) = vec_dummy(:) ! if everything is ok, take this position 
! 
!
! Now positions are replicated for the right vapour box. Translation to the second box          
!
        r0(i_part + n_chain_vap*n_mon_d/2,1) = r0(1,i_part) +  boundary(1)/2. + boundary_d(1)/2.
        r0(i_part + n_chain_vap*n_mon_d/2,2:3) = r0(2:3,i_part)  
! PBC in x and y
        if (r0(1,i_part) < 0 ) r0(1,i_part) = r0(1,i_part) + boundary(1)  ! X
        if (r0(i_part + n_chain_vap*n_mon_d/2,1) > boundary(1)  )   &
        
        r0(i_part + n_chain_vap*n_mon_d/2,1) = r0(i_part + n_chain_vap*n_mon_d/2,1) - boundary(1)  ! Y 

        if (r0(2,i_part) < 0 ) r0(1,i_part) = r0(2,i_part) + boundary(2) 
        if (r0(2,i_part) > boundary(2)  ) r0(1,i_part) = r0(2,i_part) - boundary(2) 
        if (r0(i_part + n_chain_vap*n_mon_d/2,2) > boundary(2)  )   &
        r0(i_part + n_chain_vap*n_mon_d/2,2) = r0(i_part + n_chain_vap*n_mon_d/2,2) - boundary(2) 
        if (r0(i_part + n_chain_vap*n_mon_d/2,2) < 0.  )   &
        r0(i_part + n_chain_vap*n_mon_d/2,2) = r0(i_part + n_chain_vap*n_mon_d/2,2) + boundary(2) 

     end do  ! particles in chains

end  do      ! chains   
!stop
      print*,"*     First conf for vapour phase generated ! "
!
!debug
if (debug) then
   call write_conf(1,r0(:,n_mon_tot-n_mon_vap:n_mon_tot),10)
   stop
end if
end subroutine gen_vapour
