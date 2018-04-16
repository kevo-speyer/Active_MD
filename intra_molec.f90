      subroutine intra_molec()
#include "control_simulation.h"     
      use commons 
      implicit none
      logical :: debug=.false.
      integer:: i,j
! Used to calculate 
      real (kind=8) :: l_0=0.8d0
      real (kind=8) :: r_vers(3)
      
! v_intra_molec = total intramolecular energy. ie brush + melt/droplet
! v_intra_molec_d = intramolecular energy of  melt/droplet
! v_intra_molec_e = intramolecular energy of  particle 4 chains 

           v_intra_molec = 0.
      
! Intra for brushes

    do i_chain = 1,n_chain !this goes to n_chain =>  brush only
        do i_mon = 1,n_mon-1
            i_part = (i_chain-1)*n_mon+i_mon
            j_part = i_part+1

            delta_r(1) = r0(1,i_part) - r0(1,j_part)
            delta_r(2) = r0(2,i_part) - r0(2,j_part)
            delta_r(3) = r0(3,i_part) - r0(3,j_part)

            delta_r(1) = delta_r(1) - boundary(1)*   int(2.*delta_r(1)*inv_boundary(1))
            delta_r(2) = delta_r(2) - boundary(2)*   int(2.*delta_r(2)*inv_boundary(2))
#if SYMMETRY == 1
            delta_r(3) = delta_r(3) - boundary(3)*   int(2.*delta_r(3)*inv_boundary(3))
#endif


            r_2 = 0.
            do i_dim = 1,n_dim
                r_2 = r_2 + delta_r(i_dim)* delta_r(i_dim)
            end do

#ifdef FENE_l0
            r_vers(:) = delta_r/sqrt(r_2)
#endif

!   -----  check whether interaction takes place

#ifndef NO_WARNS
            if (r_2.gt.0.98*r_chain_2)                                      &
                write(*,*) "bond stretched",i_time,r_2*inv_r_chain_2
#endif     
            r_dummy = min(r_2*inv_r_chain_2,0.98)
            pot_loc = -0.5*log(1.-r_dummy)
            v_intra_molec = v_intra_molec + pot_loc
            r_dummy = k_chain/(1.-r_dummy)*inv_r_chain_2

#ifdef FENE_l0
            delta_r(:) = delta_r(:) - l_0*r_vers(:)
#endif     

            do i_dim = 1,n_dim
                force_loc(i_dim) = r_dummy*delta_r(i_dim)
            end do
            ! HPC
            force(1,i_part) = force(1,i_part) - force_loc(1)
            force(2,i_part) = force(2,i_part) - force_loc(2)
            force(3,i_part) = force(3,i_part) - force_loc(3)

            force(1,j_part) = force(1,j_part) + force_loc(1)
            force(2,j_part) = force(2,j_part) + force_loc(2)
            force(3,j_part) = force(3,j_part) + force_loc(3)

#   if SYMMETRY == 1            
! potential part of press tensor (see viscosity.f90) ~ virial contribution
        do i = 1,3
            do j = 1,3
                press_tensor(i,j) =  press_tensor(i,j) + force_loc(i)*delta_r(j)
            end do 
        end do
# endif
        end do
    end do
    !

    v_intra_molec = k_chain*v_intra_molec

    !
    ! -----  Intramolecular potential for the melt/droplet  ****
    !
    v_intra_molec_d = 0.

    ! Warning: as it is uses the same r_chain, k_chain and interaction params,  than the brushes

    do i_chain = 1,n_chain_d !this goes to n_chaini_d =>  drop/melt only
        do i_mon = 1,n_mon_d-1
            i_part = part_init_d  + (i_chain-1)*n_mon_d+i_mon ! starts to count from the first monomer
            j_part = i_part+1
            delta_r(:) = r0(:,i_part) - r0(:,j_part)

!------ get boundaries right. PBC in x-y

#if SYMMETRY == 0
            do i_dim = 1,n_dim-1
#elif SYMMETRY == 1
!------ get boundaries right. PBC in x, y, z
                do i_dim = 1,n_dim
#endif

                    delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)*inv_boundary(i_dim))
                end do

                r_2 = 0.
                do i_dim = 1,n_dim
                    r_2 = r_2 + delta_r(i_dim)*delta_r(i_dim)
                end do

#ifndef NO_WARNS
                if (r_2.gt.0.98*r_chain_2)   write(*,*) "[intra_molec] bond stretched in the melt/droplet. t,",i_time,r_2/r_chain_2
#endif        
                r_dummy = min(r_2*inv_r_chain_2,0.98)
                pot_loc = -0.5*log(1.-r_dummy) !Fene potential ! cla ? 
                v_intra_molec_d = v_intra_molec_d + pot_loc
                r_dummy = k_chain/(1.-r_dummy)*inv_r_chain_2

                force_loc(:) = r_dummy*delta_r(:)
                force(:,i_part) = force(:,i_part) - force_loc(:)
                force(:,j_part) = force(:,j_part) + force_loc(:)

#   if SYMMETRY == 1            
! potential part of press tensor (see viscosity.f90) ~ virial contribution
        do i = 1,3
            do j = 1,3
                press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)
            end do 
        end do
# endif
            end do
        end do

        ! Intramolecular potential for the drop/melt

        v_intra_molec_d = k_chain*v_intra_molec_d


#ifdef PARTICLE_4
!
! -----  Intramolecular potential for the Particle 4 molecules
!
        v_intra_molec_e = 0.

        ! Warning: as it is uses the same r_chain, k_chain and interaction params,  than the brushes

        do i_chain = 1,n_chain_e !this goes to n_chaini_d =>  drop/melt only
            do i_mon = 1,n_mon_e-1

                i_part = part_init_e + (i_chain-1)*n_mon_e+i_mon 
                j_part = i_part+1
                delta_r(:) = r0(:,i_part) - r0(:,j_part)

! ----- get boundaries right. PBC in x-y

#if SYMMETRY == 0
                do i_dim = 1,n_dim-1
#elif SYMMETRY == 1
                    do i_dim = 1,n_dim
#endif
                        delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)*inv_boundary(i_dim))
                    end do

                    r_2 = 0.
                    do i_dim = 1,n_dim
                        r_2 = r_2 + delta_r(i_dim)*delta_r(i_dim)
                    end do
#ifndef NO_WARNS
                    if (r_2.gt.0.98*r_chain_2)    &
                        write(*,*) "[intra_molec] bond stretched in part. 4 chains",i_time,r_2/r_chain_2
#endif        
                    r_dummy = min(r_2*inv_r_chain_2,0.98)
                    pot_loc = -0.5*log(1.-r_dummy) !Fene potential 
                    v_intra_molec_e = v_intra_molec_e + pot_loc
                    r_dummy = k_chain/(1.-r_dummy)*inv_r_chain_2

                    force_loc(:) = r_dummy*delta_r(:)
                    force(:,i_part) = force(:,i_part) - force_loc(:)
                    force(:,j_part) = force(:,j_part) + force_loc(:)

#   if SYMMETRY == 1            
! potential part of press tensor (see viscosity.f90) ~ virial contribution
        do i = 1,3
            do j = 1,3
                press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)
            end do 
        end do
# endif
                end do
            end do

! Intramolecular potential for particle 4

        v_intra_molec_e = k_chain*v_intra_molec_e

#endif PARTICLE_4


#ifdef STARS

        i_part=part_init_star

        d_part=0

        do i_star = 1,n_stars !the stars

            i_part=i_part+1
            j_part=j_part+1

            do i_arm = 1,n_arms
                do i_mon_arm=1,n_mon_arm



                    if(d_part.gt.0.and.mod(d_part,n_mon_arm).eq.0) then
                        i_part=i_part-n_mon_arm*(i_arm-1)
                        j_part=j_part+1       



                    else
                        j_part = i_part+1



                    endif

                    d_part=d_part+1


                    delta_r(:) = r0(:,i_part) - r0(:,j_part)
                    !------ get boundaries right. PBC in x-y
                    do i_dim = 1,n_dim-1
                        delta_r(i_dim) = delta_r(i_dim) - boundary(i_dim)*int(2*delta_r(i_dim)/boundary(i_dim))
                    end do
                    r_2 = 0.
                    do i_dim = 1,n_dim
                        r_2 = r_2 + delta_r(i_dim)**2
                    end do
                    !-----  check whether interaction takes place
#ifndef NO_WARNS
                    if (r_2.gt.0.98*r_chain_2)   write(*,*) "bond stretched in star",i_time,r_2/r_chain_2
#endif        
                    r_dummy = min(r_2/r_chain_2,0.98)
                    pot_loc = -0.5*log(1.-r_dummy)  
                    v_intra_molec_d = v_intra_molec_d + pot_loc
                    r_dummy = k_chain/(1.-r_dummy)/r_chain_2
                    !         num_lay_fen=int((r0(3,i_part)+r0(3,j_part))/2/dz_fen)+1
                    do i_dim=1, n_dim   
                        force_loc(:) = r_dummy*delta_r(:)
                        force(i_dim,i_part) = force(i_dim,i_part) - force_loc(i_dim)
                        force(i_dim,j_part) = force(i_dim,j_part) + force_loc(i_dim)
                    end do 
                    i_part=j_part

                end do
            end do
        end do
#endif

#if PINNED==2
force(part_init_d+1, :)=force(:,part_init_d+1)-100*(r0(:,part_init_d+1)-r_start(1,:)) 
    !! applied spring force on the pinned balls
force(part_init_e+1, :)=force(:,part_init_e+1)-100*(r0(:,part_init_e+1)-r_start(2,:))
#endif


! Total intramolecular potential
       if(debug) then
               print '(a,3(2x,f15.9))',"[intra_molec]",v_intra_molec,v_intra_molec_d,v_intra_molec_e  
       end if

        v_intra_molec = v_intra_molec + v_intra_molec_d + v_intra_molec_e

end subroutine intra_molec
