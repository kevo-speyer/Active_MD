subroutine init_system()
! *  reads in system_input file the system  size and charcteristics 
! *  DYNAMIC Allocation of varibles
    use commons
#ifdef _OPENMP
      use omp_lib
#endif
 
    implicit none
#include 'control_simulation.h'    
    integer :: n_read ! for data of system_input diagnosis
    logical :: file_exists

! --- Read system_input 

    inquire(file="system_input",exist=file_exists)

    if(file_exists) then
        write(6,fmt='(/a)',advance='no') ,"   *  Reading system data from system_input file ..."
        open(file="system_input",unit=14,status="old")
        read (14,*) ;read (14,*) ;read (14,*) ;   read (14,*) ! dummy reading  

        read(14,*) n_mon,n_chain,n_cell_w_x,n_cell_w_y !  
        read(14,*) ! dummy reading
#if SYSTEM != 1 
        read(14,*,iostat=n_read) n_mon_d,n_chain_d,drop_cell_x,drop_cell_y,z_skin
#else /*SYSTEM=1: droplets */
        read(14,*,iostat=n_read) n_mon_d,n_chain_d,drop_cell_x,drop_cell_y,z_skin,z_min_d,z_max_d ! droplet stuff 
#endif
!NOTE: drop_cell_[xy] ~ droplet size
!      z_skin ~ minimum distance between bottom and wall in z coodinate in which the melt will be generated

        if(n_read/=0) print*," * Some parameter is missing in system_input file ! "
#ifdef  PARTICLE_4      
        read(14,*) !dummy
        read(14,*,iostat=n_read) n_mon_e,n_chain_e,inter4
        if(n_read/=0) then 
            print *," WARN!! parameter missing in system_input. Line for particle 4! "
           stop 
          endif
#endif
! Chemical incompatibility, parameters for non-additive potentials 
#if SOLVENT == 2 || SOLVENT == 3 
        read(14,*) !dummy
#      ifndef PARTICLE_4        
        read(14,*,iostat=n_read) delta_sig(1)
#      else
        read(14,*,iostat=n_read) delta_sig(1:3) ! three paramters if there are P4's
#      endif

        if(n_read/=0) then 
            print*," *WARN!! Parameter missing in system_input. Line: non-additive potential "
            stop
        endif
#endif 

! Wall parameters

#if WALL==2 || WALL == 3
#    ifdef  PARTICLE_4
        read(14,*) ! dummy reading
           read(14,*,iostat=n_read) a_w,sigma_w,a_w4,sigma_w4  ! Wall interaction parameters
        if(n_read/=0) then
            print*," * Parameter missing in system_input! Line: a_w, sigma_w,a_w4,sigma_w4 "
            stop 
        endif
#    else
        read(14,*) ! dummy reading
        read(14,*,iostat=n_read) a_w,sigma_w ! droplet stuff 
        if(n_read/=0) then
            print*," * Parameter missing in system_input! Line: a_w, sigma_w "
            stop 
        endif
#    endif        
#endif

#       if SYSTEM == 2  || SYSTEM == 3
               read(14,*) ! dummy reading
!              Charges for each particle type 
               read(14,*) charge(1:4)
#              ifdef STARS
                 charge(5)=0.        !if one want to read charges for stars
                                     !!STARS flag
                                     !in system_sinput, one has to change it here! attention!!!!
                 charge(6)=0.0
#              endif
               
!              Ewald sum parameters
               read(14,*) ! dummy reading
               read(14,*,iostat=n_read) alpha_c,n_kmax ! droplet stuff 
               if(n_read/=0) print*," * Some parameter is missing in system_input file ! "
#       endif

#       if SYSTEM == 3 
#              if WALL == 2 
                if (charge(1) /= 0. ) then
                print '(/a/) ', "  *  For  charged brushes with implicit walls, brush' s heads MUST have q=0"  
                stop 
                end if 
#             endif
#       endif
#ifdef STARS

        read(14,*) !dummy
        read(14,*,iostat=n_read) n_mon_arm, n_arms, n_stars
        if(n_read/=0) then 
            print *," WARN!! parameter missing in system_input. STARS! "
           stop 
          endif
#endif

#if PINNED==1||PINNED==2

        read(14,*) !dummy
        read(14,*,iostat=n_read) dist_pinned
        if(n_read/=0) then 
            print *," WARN!! parameter missing in system_input. PINNING! "
           stop 
          endif
#endif

#       ifdef POISEUILLE
! CONST_FORCE: External constant force parameters  [POISEUILLE VERSION ]
!
        read(14,*) ! dummy reading
        read(14,*,iostat=n_read) const_force(:) ! droplet stuff 
        if(n_read/=0) then
            print*," * Some parameter is missing in system_input file ! "
            stop
        end if
#       endif
        print *,"   Done!"
        close(unit=14)
    else   ! IF system_input does not exist
        print*,"'system_input' file doesn't exist in directory... "
        print*,"Finishing here."
        stop
    endif

! Write out some read parameters
#   if SYMMETRY != 1 
#   if WALL==2
        print '(/a)', "   ------- Wall fluid parameters for implicit wall--------" 
        print '(/a,2f8.4)', "  * LJ wall params: sigma_w,a_w =",sigma_w,A_w
#   ifdef PARTICLE_4        
        print '(a,2f8.4/)', "  * LJ wall params (particle 4): sigma_w,a_w =",sigma_w4,A_w4
#   endif
#   elif WALL==1
        print'(/a)', "   ------- Explicit wall--------"
#   endif
#endif /* symmetry */

#ifdef STARS
        
        print '(/a)', "   ------- STAR Parameters --------" 
        
        print '(a,3i4/)', "  * n_mon_arm, n_arms, n_stars =",n_mon_arm, n_arms, n_stars
#endif
#if SYSTEM == 2 || SYSTEM == 3 
        print '(/a)', "   ------- Coulomb Parameters --------" 
        print '(a,7f7.3)', "    Charges for each type of particle= ",charge(1:n_type) 
        print '(a,f7.3,i6/)', "    Ewald parameters: alpha, n_kmax =  ", alpha_c,n_kmax
#endif
#if SYSTEM == 1 /* droplets */
        print '(/a)', "   ------- Droplets Parameters --------" 
        print '(a,2f10.5)', "    Drop initial size and shape: drop_cell_x,drop_cell_y=   ",drop_cell_x,drop_cell_y
        print '(a,f10.5,x,f10.5)', "    Drop min and max z in the very first conf: z_min_d,z_max_d= ",z_min_d,z_max_d
#endif
!
! ----  Setting  Total number of everything    
!
       n_mon_tot =n_mon*n_chain+n_mon_d*n_chain_d  ! All the particles, but the wall.
#ifdef PARTICLE_4      
       n_mon_tot =n_mon_tot + n_mon_e*n_chain_e  ! All the particles, including type 4
       n_tot_e = n_mon_e*n_chain_e
#      ifdef STARS
       n_mon_tot = n_mon_tot + n_stars + n_stars*n_mon_arm*n_arms
#      endif
#endif

#   if WALL == 1 /* explicit wall */ 
           n_wall=4*n_cell_w_x*n_cell_w_y             ! The particles of the wall
#   else
           n_wall=0
#   endif

       n_part=  n_mon_tot + n_wall      ! All the particles of the system (wall + fluid)
       
!   Check if the algorithm combination makes sense with what it is really implemmented
       

!
!******************** Dynamic allocation of variables *********
!

! Position and forces:

       allocate (   r0(n_dim,n_part),    &
                r0_old(n_dim,n_part),    &   
                spabs(n_dim,n_part),     &       
                force(n_dim,n_part)  ) 
! Velocities and accelerations: 

        allocate ( v(3,n_part),a(3,n_part),v_half(3,n_part) ) 

#ifdef STARS
      allocate (r0_star(n_stars, n_dim))
#endif


#if STORE == 1   /* if I/O is done with unfolded coordinates  */
      allocate (r0_unfold(n_dim,n_part))
#endif
#ifdef DIFF    /* diffusion calculation */
      allocate (r0_unfold(n_dim,n_part))
#endif
                
! Reference  configuration for the diffusion coeficient
#   ifdef DIFF
       allocate ( r0_ref(n_dim,n_part) )
#   endif       
! Vector used for the thetered chains heads (used if not explicit wall)               
!obs       allocate( r_head_old(n_chain*n_mon,3) ) !alpha version
                
! Equilibrium sites of wall atoms:
       allocate(r_wall_equi(n_dim,n_wall))

! veloc., accel. and higher derivatives
!obs       allocate(rx(n_part,n_dim,n_order))  
! Mass and atom types, vectors:
        allocate(   a_type(n_part),              &
                     mass(n_part),               &
                     inv_mass(n_part),               &
                     random(n_part)   )         ! random vector
! Variables to monitor minimum image convention 
        allocate  ( mic_count(n_dim,n_part) )
! Variables used for binning:

                   n_bin_x= int(0.4*n_cell_w_x)
                   n_bin_y= int(0.7*n_cell_w_y)
                   n_bin_z= 20   !ori: 20                                     

! NOTE: n_neigh is the maximum number of neighbors of each atom
! it is defined to a more or less arbitrary value of 256 or 300
! It is also a typical reason for which the program crashes. 

allocate (  ff_list(0:n_neigh_fl+n_layer*n_wall,n_part) ) ! ,             &
!              fw_list(0:n_neigh_wa,n_mon_tot),                          &
!              ww_list(0:n_neigh_ww,n_wall/2)      )
!print *,size(ff_list,dim=2) ; stop
       allocate( mic_old(n_dim,n_part) )



! Variables for tensor of gyration
!obs        allocate (    chcm(1:n_chain,1:n_chain),                            &
!obs                    gyrten(1:n_dim,1:n_dim,1:10),                         &
!obs                    gyrt(1:n_dim,1:n_dim,1:n_chain)  )  

!cla: force on brush heads variable and CM diffusion vector
        allocate ( f_on_heads(3,n_chain) ) !, r0_cm_m(n_chain_d,3) 

! Variables for the DPD implemmentation

        allocate (force_d(3,n_mon_tot) , force_r(3,n_mon_tot) )

!TO OPRIMIZE PARALLEL VERSION
#ifdef _OPENMP
!$OMP PARALLEL 
        numth=omp_get_num_threads() !Define parellel variables for each thread
!$OMP END PARALLEL
!        allocate (force_LJ(3,n_mon_tot,0:numth-1))
!        allocate (v_fl_fl_tmp(0:numth-1))
#endif
!/TO OPRIMIZE PARALLEL VERSION



!        allocate(random_2(n_part),random_3(n_part))

! Variables for coulomb calculations 

#if SYSTEM == 2  || SYSTEM == 3  
        allocate ( q(n_mon_tot) ) 
#endif

end subroutine init_system
