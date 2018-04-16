     subroutine init_params

! This routine:
! - sets simulation parameters if mfa_input doesn't exist
! - reads simulation parameters from mfa_input if it exists
! -  defines cutoffs and shifts in energy
! - writes out the simulation parameters for a new run and aditional info
!   in mfa_output
! - initialize predictor/corrector or verlet integration algorhitms
! - Initializes other parameters as well.
      use commons

#ifdef _OPENMP
      use omp_lib
      use Par_Zig_mod 
#endif

      use ziggurat ! gaussian number generator suite
#include 'control_simulation.h'
      implicit none
      real (kind=8) :: init_dens,graft_dens,wall_t_area
      integer :: iseed_z,j 
      logical :: input_file,ref_file,es
      character (len=6) :: dummy_char
!   
!---  read parameters that determine experiment
        inquire(file="mfa_input",exist=input_file)

        if(.not.input_file) then    ! DEFAULT PARAMETERS  

            write(*,*) "not adapted for default parameters! Check!!"
            stop
            print '(/a/)'," *     'mfa_input' file not found ... , using default parameters"
            s_time = 1
            n_relax = 50000     !ori: 1000
            n_obser = 50000     !ori: 1000
            n_safe  = 1000         !ori: 2000
            dt = 0.0001 !ori: 0.005
            !obsolete       temp_time =  0.
            !obsolete       temp_init = 1.68 !1.68 !wanted=1.68  !ori: 0.5
            !obsolete       temp_final = 1.68 !1.68 !ori: 0.5

            temp = 1.68 ! default temp definition

            !obsolete f_wall_fix = 1 !! que es ?
            x_space = 1.2091356
            !      force pressure conv factor: 1.2557
            x_space = 1.05
            i_cell_w_x = n_cell_w_x
            i_cell_w_y = n_cell_w_y
            if((i_cell_w_x.eq.i_cell_w_y).or.                                &
                (i_cell_w_x.eq.2*i_cell_w_y)) then
!----   commensurate surfaces
            y_space = x_space*1.7320508                     ! this factor cla?
        else if(((i_cell_w_x.eq.31).and.(i_cell_w_y.eq.18)).or.          &
            ((i_cell_w_x.eq.62).and.(i_cell_w_y.eq.36))) then
!----   incommensurate surfaces
        y_space = 31.*x_space/18.                         ! cla?
    else
        write(*,*) "look up dimensions of walls"
        stop
    end if   ! finishes default parameters list
!       z_space_wall = 2.62*n_mon_tot/n_wall+2.5
    z_space_wall = 28. ! !ori=50 !Torsten: typical value 17.5 ! for n=10 chains = 30.

! Default mass definitions 

    mass_type(1) = 1.
    mass_type(2) = 1.
    mass_type(3) = 1.
!    if(f_explicit_wall) then ! IF IT IS THE ORIGINAL PROGRAM WITH explicit walls
#       if WALL == 1
        mass_type(4) = 15.
#    else /* not explicit walls */
        mass_type(4) = 1.
#    endif

    epsil(1,1) = 1.
    epsil(1,2) = 1.
    epsil(1,3) = 1.
!obsol       if(ad_flag.eq.0) then    ! => adsorbed, not grafted
!obsol           epsil(1,4) = 1.
!obsol       else
    epsil(1,4) = 100. 
!obsol       end if
    epsil(2,1) = 1.
    epsil(2,2) = 1.
    epsil(2,3) = 1.
    epsil(2,4) = 1.
    epsil(3,1) = 1.
    epsil(3,2) = 1.
    epsil(3,3) = 1.
    epsil(3,4) = 1.
    !obso       if(ad_flag.eq.0) then   ! adsorbed, not grafted
    !obso           epsil(4,1) = 1.
    !obso       else
    epsil(4,1) =100. 
    !obs       end if
    epsil(4,2) = 1.
    epsil(4,3) = 1.
    epsil(4,4) = 1.
    sigma(1,1) = 1.
    sigma(2,2) = 1.
    sigma(3,3) = 1.
    sigma(4,4) = 1.
    !Just in case, confirm values of epsilon if no explicit walls
!    if(.not.f_explicit_wall) then
#       if WALL == 2 || WALL == 3 || WALL == 4        
        epsil(:,:) = 1.
        !There is not interaction between brush heads
        epsil(1,1) = 0.
#       endif        
!    end if

    friction(:) = 0.5
    k_spring_wall = 100.             ! cla?
    !obsolete f_cut_off = 0                    ! for first conf high dens whether uses or not cutoff 
    iseed = 311499519
    k_chain = 30.
    r_chain = 1.5

    ! ---- Switching on of the force 


    f_twall(1) = 2 ! behavior of the top wall. xyz directions cla?
    f_twall(2) = 2   ! so ft_wall(:) = 1 => v=0 for the 3 coor. 
    f_twall(3) = 1
    va_spring_twall(1) = 0.
    va_spring_twall(2) = 0.
    va_spring_twall(3) = 0.
    k_spring_twall(1) = 0.
    k_spring_twall(2) = 0.
    k_spring_twall(3) = 0.
    !       f_minimize = 0               !cla?
    !       n_list_up = 1                !cla?

    !------------- end default params -------------------------------------------       

else      ! ---- READ mfa_input 

    open(unit=10,file="mfa_input",status="old")
    write (6,fmt='(/a/)',advance='no') "  *   Reading parameters from mfa_input ..."
    read(10,*)   s_time     ! I don't know if it will be used 
    read(10,*) n_relax         ;print '(/a30,i6)',  "n_relax ",n_relax     
    read(10,*) n_obser         ;print '(a30,i6)',   "n_obser ",n_obser        
    read(10,*) n_safe          ;print '(a30,i6)',   "n_safe  ",n_safe             
    read(10,*) dt              ;print '(a30,g12.5)',"dt      ",dt                
    read(10,*) !c_dummy        
    read(10,*) !r_dummy       
    read(10,*) temp          ;  print '(a30,g16.5)',"Temperature  ",temp                  
    read(10,*) 
    read(10,*) 
    read(10,*) 
    read(10,*) x_space         ;  print '(a30,g12.5)' ,"   x_space      ",x_space
    read(10,*) y_space         ;  print '(a30,g12.5)' ,"   y_space      ",y_space               
    read(10,*) z_space_wall    ;  print '(a30,g12.5)' ,"   z_space_wall ",z_space_wall          
    read(10,*) !c_dummy         !;print*,"c_dummy     ",c_dummy           

! ---- Read interaction parameters 

    do i_type = 1,n_type  ! Reading loop 
        read(10,*) mass_type(i_type) ; print'(a30,g12.5)',      "mass_type =", mass_type(i_type)
        do j_type=1,n_type 
            read(10,*) epsil(i_type,j_type) ; print '(a30,g12.5)',  "epsil(i,j) =", epsil(i_type,j_type)
        end do
        read(10,* ) sigma(i_type,i_type)  ; print '(a30,g12.5)',"sigma(i,i) =", sigma(i_type,i_type) 
        read(10,*) friction(i_type),c_dummy ; print '(a30,g12.5)',"friction =" ,friction(i_type)
    end do

! ---- Redifinition of interaction parameters [[check]]  

        do j_type = 1,n_type
#   if WALL == 1   
#       ifdef PARTICLE_4
                epsil(n_type,:)=1.
#          ifdef FREE_HEADS
                epsil(1,n_type)=250.
                epsil(n_type, 1)=250.
#          endif
                cycle
            endif
#     endif /* particle 4*/

            if (j_type.eq.n_type) then
                sigma(n_type, n_type)=1.
                mass_type(n_type)=1. ! WARN: ori 150
                friction(n_type)=friction(n_type-1)
                cycle
            endif

#   endif /*wall=1: explicit wall*/              

#       ifdef STARS
            if ((j_type.eq.5).or.(j_type.eq.6)) then
                epsil(:,5)=1. !! the interaction with stars is fixed. Epsilon is standard for STARS
                epsil(5,:)=1.
                epsil(:,6)=1. 
                epsil(6,:)=1.

                sigma(5,5)=1.
                sigma(6,6)=1.
                mass_type(5)=1.

                mass_type(6)=1.
                friction(5)=friction(4)
                friction(6)=friction(4)
                cycle
            endif
#       endif

        end do
     

    read(10,*) k_spring_wall,c_dummy ; print '(/a30,g12.5)', "   k_spring_wall = ", k_spring_wall
    
# ifdef BENDING        
read(10,*) k_bend ; print '(a30,e14.6)',"k_bend = ", k_bend !bending elastic constant, set default to 0
read(10,*) alpha_eq ; print '(a30,e14.6)',"alpha_eq = ", alpha_eq  ! equilibrium bending angle set default to 0
! Active matter parameters. It must be beta1<beta2, no checking logic for now.
read(10,*) beta1 ; print '(a30,e14.6)',"beta1 (deg) = ", beta1  ! Angle measured from horizontal plane
read(10,*) beta2 ; print '(a30,e14.6)',"beta2 (deg) = ", beta2  ! Angle measured from horizontal plane.
# endif

read(10,*) k_act ; print '(a30,e14.6)',"k_act = ", k_act  ! Active k constant used in active brush (metronome)
read(10,*) k_spr_x ; print '(a30,e14.6)',"k_spr_x = ", k_spr_x  ! Spring constant, x direction (spring_array)
read(10,*) k_spr_y ; print '(a30,e14.6)',"k_spr_y = ", k_spr_y  ! Spring constant, y direction (spring_array)
    
    read(10,*) !i_dummy
    read(10,*) !c_dummy
    read(10,*) iseed ; print '(a30,i10)',"iseed = ",iseed
    read(10,*) !c_dummy
    read(10,*) k_chain ; print '(/a30,f12.7)',"  FENE: kchain = ",k_chain
    read(10,*) r_chain ; print '(a30,f12.7)',"  FENE: rchain = ",r_chain
    read(10,*)! c_dummy
    read(10,*)! c_dummy
    read(10,*)! c_dummy
    read(10,*)! c_dummy
    read(10,*) (f_twall(i_dim),i_dim=1,n_dim) ; print '(a30,3i6)', "    f_twall =",f_twall(:)
    read(10,*) (va_spring_twall(i_dim),i_dim=1,n_dim) ; print '(a30,3f10.5)', "    Wall velocity = ",va_spring_twall(:)
#ifndef SHEARED
    read(10,*) (k_spring_twall(i_dim),i_dim=1,n_dim) ;print '(a30,3f10.5)', "    k_spring_twall(:) = ",k_spring_twall(:)
#else
    read(10,*) (turn_time(i_dim),i_dim=1,n_dim) ; print '(a30,3f10.5)',"    turn_time=" , turn_time(:)
    k_spring_twall(:) = 0.0
    read(10, *) p_time
#endif
!    read(10,100) c_dummy
!    read(10,110) ! obsolete f_minimize
!    read(10,*) ! n_list_up
    print '(/a/)',"           Done !"
    close(10)
end if

!---  write out mfa_output
!
       inquire (file='mfa_output',exist=ref_file)
       if (ref_file) then
           print '(/a/)' ,"  *  WARNING: Overwriting mfa_output file " 
       end if
       open(20,file="mfa_output",status="unknown")

       write(20,210) s_time+n_relax+n_obser , " initial integer starting time"
       write(20,210) n_relax," number of relaxation steps"
       write(20,210) n_obser," number of observation steps"
       write(20,210) n_safe ," # of steps between configuration storage"
       write(20,201) dt     ," time step increment"
       write(20,*)
       write(20,'(f12.5,a)') temp ," obsolete"
       write(20,'(f12.5,a)') temp ," temperature"
       write(20,'(f12.5,a)') temp ," obsolete"
       write(20,*)
       write(20,*) "  obsolete: flagg for constraining wall atoms"
       write(20,201) x_space      ," spacing betw wall units in x direc"
       write(20,201) y_space      ," spacing betw wall units in y direc"
       write(20,201) z_space_wall ," initial inter-wall spacing"
       write(20,*)

       do i_type = 1,n_type
           
#if WALL==1   
#    ifdef PARTICLE_4
      if (i_type.eq.n_type) cycle
#    endif
#endif
#ifdef STARS

           if ((i_type.eq.5).or.(i_type.eq.6)) cycle
#endif
           write(20,201) mass_type(i_type)   ," mass"
           do j_type = 1,n_type

#if WALL==1   
#    ifdef PARTICLE_4
      if (j_type.eq.n_type) cycle
#    endif
#endif

#ifdef STARS

           if ((j_type.eq.5).or.(j_type.eq.6)) cycle
#endif
               write(20,201) epsil(i_type,j_type)," Lennard Jones epsilon"
           end do
           write(20,201) sigma(i_type,i_type)," Lennard Jones sigma"
           write(20,201) friction(i_type)    ," friction constant"
       end do
       write(20,201) k_spring_wall," spring const. to wall equil. sites"
       
# ifdef BENDING        
       write(20,201) k_bend," bending elastic constant"  !bending elastic constant, set default to 0
       write(20,201) alpha_eq," bending equilibrium angle"  ! equilibrium bending angle set default to 0
       write(20,201) beta1," beta1 (deg), active parameter"  ! Angle measured from horizontal plane
       write(20,201) beta2," beta2 (deg), active parameter"  ! Angle measured from horizontal plane.
# endif

       write(20,201) k_act," k_act, elastic constant, active parameter"  ! 
       write(20,201) k_spr_x," k_spr_x, elastic constant, spring_array x direction"  ! 
       write(20,201) k_spr_y," k_spr_y, elastic constant, spring_array x direction"  ! 

       write(20,210) !f_cut_off," cut-off flag (0) short (1) long range [OBSOLETE]"
       write(20,*)
       write(20,210) iseed+2," seed for random number generator"
       write(20,*)
       write(20,201) k_chain," effective spring constant within chain"
       write(20,201) r_chain," max. spacing between neighbored monomers"
       write(20,*)
       write(20,*)
       write(20,'(a)') "  flaggs and variables controlling boundary conditions"
       write(20,*)
       write(20,130) (f_twall(i_dim),i_dim=1,n_dim)
       write(20,103) (va_spring_twall(i_dim),i_dim=1,n_dim)
#ifndef SHEARED
       write(20,103) (k_spring_twall(i_dim),i_dim=1,n_dim)
#else

       write(20,130) (turn_time(i_dim),i_dim=1,n_dim)
       write(20,110) p_time
#endif
       write(20,*)
       write(20,*) " obsolete: f_minimize flagg to relax to next minimum"
       write(20,*) " obsolete: n_list_up ,  # of steps between list up date"
!
!---  Extra information about parameter statements
!
      write(20,*)
      write(20,*) "----------  Additional info. Not used as input  ----------------"
      write(20,*)
      write(20,210) n_mon     ,"  monomers per chain"
      write(20,210) n_chain   ,"  number of 'brush' chains"
      write(20,210) n_layer*n_wall   ,"  number of fluid particles"
      write(20,210) n_chain_d  ,"  number of droplet/melt chains"
      write(20,210) n_chain_d*n_mon_d   ,"  number of droplet/melt particles"
      write(20,220) n_cell_w_x,n_cell_w_y,"  number of wall cells in x-y direction"
      write(20,*)
      write(20,210) n_order   ,"  order of integrator"
      write(20,*) "    obsolete: ad_flag  flag for adsorbed or endgrafted"
!      if(ad_flag.eq.0) then
!          write(20,210) ad_flag   ,"    adsorbed polymers"
!      end if
!!      if(ad_flag.eq.1) then
          write(20,*) "      endgrafted polymers"
!!      end if
      write(20,*) "   "
!      if(f_explicit_wall) then
#       if WALL == 1 
          write(20,'(a)') " * Using explicit wall atoms with LJ int."
#      else
          write(20,'(a)') " * Using potential 9-3, for wall int."
#      endif
      ! informs if DPD or not      
      write(20,*) "   "

!      if(f_dpd) then
# IF THERMOSTAT == 0 
          write(20,'(a/)') " * Using DPD thermostat and Veloc. Verlet integration scheme."
          write(20,'(a,f8.5)') " *       DPD Friction value in all directions: gamma = ",friction(1)
          print '(a,f8.5)', " *       DPD Friction value in all directions: gamma = ",friction(1)
          print '(a/)'    , "        ( Taken from the first friction value in mfa_input )   "        
!      else
# elif THERMOSTAT == 1
          write(20,'(a50)') " * Using Langevin thermostat and predictor corrector integrator  "  
          write(20,'(a50)') "   (If integrator order = 2 ~ velocity verlet)  " 
# endif
!      end if
      write(20,*) "   "
      !
!obsolete      if((temp_time.ge.n_relax*dt).and.(temp_init.ne.temp_final))       &
!obsolete      write(*,'(a)') "Temperature ramp might take too long."
!
! ---- Write out  more info
!
        print '(/a,4f10.5/)', "mass_type= ",mass_type(:)

      init_dens = dble(n_mon_tot)/(z_space_wall*dble(n_cell_w_x*n_cell_w_y)*x_space*y_space)

   !grafting density [1/sigma^2]
      wall_t_area = dble(n_cell_w_x*n_cell_w_y)*x_space*y_space ! assume this is en units of sigma
      graft_dens = dble(n_chain)/(wall_t_area*2.0) ! divide the numbers of chains over the  
      ! total wall area
      write(20,"(3f10.5,a50)") n_cell_w_x*x_space,n_cell_w_y*y_space,z_space_wall, " MD Box dimensions"  
      write(20,'(g18.9,a50)') init_dens, "      initial monomer density "
      write(20,'(2f18.5,a50)') wall_t_area,graft_dens, "one wall area, grafting density     "
!      
!---  Initialize variables that immediatly follow from input variables
!
      dt_2 = dt**2
      r_chain_2 = r_chain**2
      inv_r_chain_2 = 1./r_chain_2
      k_chain = k_chain*r_chain_2
!obso      if(f_wall_fix.eq.1) then
!obso          n_loop = n_mon_tot
!obso      else
!obso          n_loop = n_part
!obso      end if

!---  Apply standard sum rules to Lennard Jones parameter
! Additive potential 
#if SOLVENT == 0 || SOLVENT == 1
      do i_type = 1,n_type-1
          do j_type = i_type+1,n_type
              sigma(i_type,j_type) = (sigma(i_type,i_type)+sigma(j_type,j_type))/2
              sigma(j_type,i_type) = sigma(i_type,j_type)
          end do
      end do
#       ifndef PARTICLE_4      
#       if WALL != 1 
              sigma(:,n_type) = 0.
              sigma(n_type,:) = 0.
#       endif
#       endif      
#endif      

#if SOLVENT == 2 || SOLVENT == 3

! NON-additive potential
! First apply standard sum rules and then correct with the input data
      do i_type = 1,n_type-1
          do j_type = i_type+1,n_type
              sigma(i_type,j_type) = (sigma(i_type,i_type)+sigma(j_type,j_type))/2
              sigma(j_type,i_type) = sigma(i_type,j_type)
          end do
      end do
! brush-melt interaction is non- additive !      
! interaction brush-melt 
      sigma(2,3) = delta_sig(1)
      sigma(3,2) = delta_sig(1)
      sigma(1,3) = delta_sig(1)
      sigma(3,1) = delta_sig(1)

#       ifdef PARTICLE_4      

! interaction particle4-brush

      sigma(2,4) = delta_sig(2)
      sigma(4,2) = delta_sig(2)
      sigma(1,4) = delta_sig(2)
      sigma(4,1) = delta_sig(2)

! interaction particle4-liquid

      sigma(3,4) = delta_sig(3)
      sigma(4,3) = delta_sig(3)

#       endif

#endif      
      print '(/a/)' , "  -------- LJ Interaction parameters -------     " 
      do i_type = 1,n_type
          write(*,'(a,i4, 7f10.4)') "   Sigma: particle   ",i_type,sigma(i_type,:)
      end do
      print '(a/)'," "      
      do i_type = 1,n_type
          write (*,'(a,i4, 7f10.4)') "   Epsilon: particle ",i_type,epsil(i_type,:)
      end do
!pause
!---  Absorb factor 4 into epsilon and initialize sigma_2
!
      do i_type = 1,n_type
          do j_type = 1,n_type
              sigma_2(i_type,j_type) = sigma(i_type,j_type)**2
              epsil(i_type,j_type) = 4.*epsil(i_type,j_type)
          end do
      end do


! Additional variables to improve performance
! NOTE: When we have implicit walls with fixed heads of the brushes we don' t
! take the head for averages so they are not included in 1/N

    inv_dt = 1./dt
    inv_dt_2 = inv_dt**2
#if (SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3 ) && SYMMETRY == 0
    inv_N = 1./dble(n_mon_tot - n_chain)
#endif
#if SYSTEM == 2 || SYMMETRY == 1
    inv_N = 1./dble(n_mon_tot)
#endif

! WARN      
!note: these definitions are not clear 
! I am overriding them, at least for bad solvent and flat wall,
! after this piece of code.

#        if WALL != 1        
            epsil(1,1) = 0.
#       endif             

!---  define cutoffs and shifts in energy

                range_2(:,:) = (2.**(1./6.)*sigma(:,:))**2

!                 if(ad_flag.eq.0) then
!
!                  range_2(1,n_type) = (2.2*sigma(1,n_type))**2 
!                  range_2(n_type,1) = (2.2*sigma(n_type,1))**2

!                  range_2(2,n_type) = (2.2*sigma(2,n_type))**2
!                  range_2(n_type,2) = (2.2*sigma(n_type,2))**2
!
!
!                    end if ! ad_flag=0 no endgrafted chains


                  range_2(1,n_type) = (2.2*sigma(1,n_type))**2 
                  range_2(n_type,1) = (2.2*sigma(n_type,1))**2 


! original cutoffs and shifts for atomic walls

!    if( f_explicit_wall ) then
#   if WALL == 1  /* explicit wall */
#           if SOLVENT == 0 || SOLVENT == 3


        range_2(:,:) = (2.*2.**(1./6.)*sigma(2,2))**2  

        range_2(1,n_type) = 0.
        range_2(n_type,1) = 0.
        range_2(n_type,n_type) = 0.


#           elif SOLVENT == 1 || SOLVENT == 2
        range_2(1:n_type-1,1:n_type-1) = (2.**(1./6.) *sigma(1:n_type-1,1:n_type-1))**2    
#           endif                    

#   endif



#       if WALL == 2 || WALL == 3 || WALL == 4
!
! ---- For all the walls, but the explicit ones 
!

#               if SOLVENT == 0 /* poor solvent */

! ----  The  cut off is  2*(2)**(1/6)* sigma  for poor solvent:        

                        range_2(:,:) = ( 2.*(2.**(1./6.))*sigma(:,:) )**2 

#                       ifndef PARTICLE_4
                        range_2(:,n_type) = 0.
                        range_2(n_type,:) = 0.
#                       endif

#               endif

! Interaction ranges for non-additive situations 
! NOTE: for SOLVENT=2 everything is repulsive but sig_23 is different and will
! have different Rc. No explicit change is needed. 

#       if SOLVENT == 3

! Here sig_23 inteacts with non-additive sigma but also interaction 23 id of the
! good solvent. Everybody else keeps being poor solvent.

! Define ranges and shifts for interactions 2-3 

! first all poor solvent 
                    range_2(:,:) = ( 2.*(2.**(1./6.))*sigma(:,:) )**2 
! redefine add good solvent 


#       endif /* solvent 3: non-additive potential */

!!!!!! Hydrpphobia makes interaction between brush and melt purely repulsive!!!!!!
# ifdef HYDROPHOBIA
range_2(1,3) = (2.**(1./6.)*sigma(1,3) )**2
range_2(2,3) = (2.**(1./6.)*sigma(2,3) )**2
range_2(3,1) = (2.**(1./6.)*sigma(3,1) )**2
range_2(3,2) = (2.**(1./6.)*sigma(3,2) )**2
# endif /* hydrophobia */

# ifdef BRUSH_IN_GOOD_SV
range_2(1,2) = (2.**(1./6.)*sigma(1,2) )**2
range_2(1,1) = (2.**(1./6.)*sigma(1,1) )**2
range_2(2,2) = (2.**(1./6.)*sigma(2,2) )**2
range_2(2,1) = (2.**(1./6.)*sigma(2,1) )**2
# endif /*brush_in_good_solvent*/

#   endif /* implicit walls */

!            end if ! not f_explicit_wall

#if WALL==2 
     do i_type=1, n_type
          sigma_wall(i_type)=(sigma_w+sigma(i_type,i_type))/2
          a_wall(i_type)=a_w
     enddo
#    ifdef PARTICLE_4
          sigma_wall(4)=(sigma_w4+sigma(4,4))/2.
          a_wall(4)=a_w4

#    endif
#if SYMMETRY != 1          
        print '(/a/)','  *  Interaction with implicit walls:' 
        print '(a,4(f8.3,x))','   sigma_w= ',sigma_wall(:)
        print '(a,4(f8.3,x))','   a_w    = ',a_wall
#endif        

#endif

! Interaction for particle 4 
#ifndef PARTICLE_4

!Just in case Rcut = 0 for wall particles and the rest of the world    

! Interaction with wall particles 

#     if WALL==2 /* implicit wall */
             range_2(:,n_type) = 0. 
             range_2(n_type,:) = 0.
#     endif
!!! already done #     if WALL==1 /* explicit wall */
!!! already done #            if SOLVENT==1 
!!! already done              range_2(:,n_type) = (2.**(1./6.)*sigma(:,n_type))**2
!!! already done              range_2(n_type,:) = (2.**(1./6.)*sigma(n_type,:))**2
!!! already done #            endif
!!! already done #            if SOLVENT==0 
!!! already done              range_2(:,n_type) = (2.2*sigma(:,n_type))**2
!!! already done              range_2(n_type,:) = (2.2*sigma(n_type,:))**2
!!! already done #            endif 
!!! already done #     endif
 
#else           /*ifdef particle 4*/ 
 
#     if WALL==1 /* explicit wall */
#            if SOLVENT==1 
             range_2(:,n_type) = (2.**(1./6.)*sigma(:,n_type))**2
             range_2(n_type,:) = (2.**(1./6.)*sigma(n_type,:))**2
#            endif
#            if SOLVENT==0 
             range_2(:,n_type) = (2.2*sigma(:,n_type))**2
             range_2(n_type,:) = (2.2*sigma(n_type,:))**2
#            endif 
#     endif


        if (inter4 == 1 ) then  ! Particle 4 has repulsion with all the other particles and with itself

            range_2(1:3,4) = (2.**(1./6.)*sigma(1:3,4))**2 ! repulsive with brush and melt 
            range_2(4,1:3) = (2.**(1./6.)*sigma(4,1:3))**2 ! repulsive with brush and melt
            range_2(4,4)   = (2.**(1./6.)*sigma(4,4))**2   ! repulsive among particles 4

           elseif (inter4 == 2 ) then ! Particle 4 has atraction among themselves and repulsion with everybody else

            range_2(1:3,4) = (2.**(1./6.)*sigma(1:3,4))**2 ! repulsive with brush and melt 
            range_2(4,1:3) = (2.**(1./6.)*sigma(4,1:3))**2 ! repulsive with brush and melt
            range_2(4,4)   = ( 2.*(2.**(1./6.)*sigma(4,4)) )**2   ! atractive among particles 4

          elseif (inter4 == 3 ) then ! Particle 4 has atraction among itself.  Particle 3 also.


! Interactions for particle 4 

            range_2(1:3,4) = (2.**(1./6.)*sigma(1:3,4))**2 ! repulsive with brush and melt 
            range_2(4,1:3) = (2.**(1./6.)*sigma(4,1:3))**2 ! repulsive with brush and melt

            range_2(4,4)   = ( 2.*(2.**(1./6.)*sigma(4,4)) )**2   ! atractive among particles 4

!
! Redefine shift for particle 3-3
!
!obsolete?                r_dummy = (sigma(3,3)**2/range_2(3,3))**3
!obsolete?                e_shift(3,3) =   r_dummy*(r_dummy-1.) 

          elseif (inter4 == 4 ) then ! Particle 4 has atraction with all the others particles 
! Interactions for particle 4 
            range_2(1:3,4) = (2.*2.**(1./6.)*sigma(1:3,4))**2 ! atractive with brush and melt 
            range_2(4,1:3) = (2.*2.**(1./6.)*sigma(4,1:3))**2 ! atractive  with brush and melt
            range_2(4,4)   = ( 2.*(2.**(1./6.)*sigma(4,4)) )**2   ! atractive among particles 4


       end if ! interaction type of particle 4
!
! e_shift definition for particle 4:
!obsolete?          do i_type=1,3      
!obsolete?                r_dummy = (sigma(i_type,4)**2/range_2(i_type,4))**3
!obsolete?                e_shift(i_type,4) = r_dummy*(r_dummy-1.)
!obsolete?                e_shift(4,i_type) = r_dummy*(r_dummy-1.) 
!obsolete?            end do

                r_dummy = (sigma(4,4)**2/range_2(4,4))**3
!obsolete                e_shift(4,4) =   r_dummy*(r_dummy-1.) 


!obsolete? !! Redefine the e-shifts for walls, if particle_4
!obsolete?         do i_type = 1,n_type
!obsolete? !            do j_type = 1,n_type
!obsolete?                 r_dummy = (sigma(i_type,n_type)**2/range_2(i_type,n_type))**3
!obsolete? !                e_shift(i_type,n_type) = r_dummy*(r_dummy-1.)
!obsolete?                  
!obsolete?                 r_dummy = (sigma(n_type,i_type)**2/range_2(n_type,i_type))**3
!obsolete?                 e_shift(n_type,i_type) = r_dummy*(r_dummy-1.)
!obsolete?          end do
!obsolete? !        end do


! NOTE:  Here more interaction definitions among different species can be added for
! other values of inter4 

#endif      /* particle 4 */      

!            
! --- Compute shifts in the potential energy
! NOTE: this block MUST go AFTER the complete definition of: sigma and range_2
!

#if SYSTEM == 0 || SYSTEM == 1|| SYSTEM==3 /* if there are brushes, we kill head head interaction*/

#ifndef FREE_HEADS
            range_2(1,1) =0.
#else
range_2(1,1) = ( (2.**(1./6.))*sigma(1,1) )**2
range_2(1,n_type)=( 2.*(2.**(1./6.))*sigma(1,n_type) )**2
range_2(n_type, 1)=( 2.*(2.**(1./6.))*sigma(n_type,1) )**2
#endif
#endif
        do i_type = 1,n_type
            do j_type = 1,n_type

                if( range_2(i_type,j_type) > 0.0 ) then 
                    r_dummy = (sigma(j_type,i_type)**2/range_2(j_type,i_type))**3
                    e_shift(j_type,i_type) = r_dummy*(r_dummy-1.)
                else
                     e_shift(j_type,i_type) = 0.0
                end if

            end do
        end do


! Here all the interactions a correct cut-offs should be defined 
!ite(*,*) sqrt( range_2(:,:)), e_shift(:,:)
         print '(/a,7f10.5)' , "  * Particle 1 interacts with cutoffs: ",sqrt( range_2(:,1) ) 
         print '(a,7f10.5)'  , "  * Particle 2 interacts with cutoffs: ",sqrt( range_2(:,2) ) 
         print '(a,7f10.5)'  , "  * Particle 3 interacts with cutoffs: ",sqrt( range_2(:,3) ) 
         print '(a,7f10.5/)' , "  * Particle 4 interacts with cutoffs: ",sqrt( range_2(:,4) ) 
#   if WALL==1
#        if SYSTEM ==3 || SYSTEM==0 
#           ifdef  PARTICLE_4
         print '(a,7f10.5/)' , "  * Particle 5 interacts with cutoffs: ",sqrt( range_2(:,5) ) 
#           endif
#        endif
#   endif

#ifdef STARS
         print '(a,7f10.5)' , "  * Particle 6 interacts with cutoffs: ",sqrt( range_2(:,6) ) 
         print '(a,7f10.5)' , "  * Particle 7 interacts with cutoffs: ",sqrt( range_2(:,7) ) 
#endif
         print '(/a,7f10.5)' , "  * Particle 1 potential's shifts: ",4.*e_shift(:,1) 
         print '(a,7f10.5)'  , "  * Particle 2 potential's shifts: ",4.*e_shift(:,2)  
         print '(a,7f10.5)'  , "  * Particle 3 potential's shifts: ",4.*e_shift(:,3)  
         print '(a,7f10.5)' , "  * Particle 4 potential's shifts: ",4.*e_shift(:,4)  
#   if WALL==1
#        if SYSTEM ==3 || SYSTEM==0 
#           ifdef  PARTICLE_4
         print '(a,7f10.5/)' , "  * Particle 5 interacts with cutoffs: ",sqrt( range_2(:,5) ) 
#           endif
#        endif
#   endif

#ifdef STARS
         print '(a,7f10.5)' , "  * Particle 6 potential's shifts: ",4.*e_shift(:,6)  
         print '(a,7f10.5/)' , "  * Particle 7 potential's shifts: ",4.*e_shift(:,7)  
#endif
# ifdef HYDROPHOBIA
 print*,""
print*, "     * Interactions between brush and melt are purely repulsive   "
print*,""
# endif /* hydrophobia */

# ifdef BRUSH_IN_GOOD_SV
 print*,""
print*, "     * Interactions between grafted polymers are purely repulsive   "
print*,""
# endif /* brush in good solvent */


! ----  Max interaction range (will be used for the binning boxes) :
!
! NOTE: we allow now for the possibility of having a DPD cut-off greater than the
! conservative forces cutoff. Improves DPD efficiency ?

        r_cut_max = sqrt( maxval(range_2(:,:)) )

        print '(/a,f10.5/)', "  * Maximum cut-off radius for potentials: ",r_cut_max

#ifdef DPD_CUT_OFF
!
!       Get DPD cut-off and redefine the cut_off_max 
!
        r_cut_dpd =  DPD_CUT_OFF
        r_cut_max = max(r_cut_dpd,r_cut_max)
        r_cut_dpd_2 =  r_cut_dpd**2
        r_cut_max_2 = r_cut_max**2 

        print '(/a,f10.5/)', "  * DPD cut-off radius for potentials: ",r_cut_dpd

#endif
        print '(/a,f10.5/)', "  * Max  cut-off radius for ALL: ",r_cut_max

! --- Aux variable for HPC

    do i_type=1,n_type
         do j_type=1,n_type
         if (range_2(i_type,j_type) /= 0 ) then
            inv_range_2(i_type,j_type) = 1./ range_2(i_type,j_type)
        else
            inv_range_2(i_type,j_type) = 0.
        end if

        end do
    end do

! --- Set skin and variables for verlet lists (?) 

      mic_count(:,:) = 0
      skin = 0.4
      print *, "  * Skin for verlet lists = ",skin


! now I take the skin as that of the maximum cut off from of all the
! interactions 

      r_dummy = r_cut_max
      skin_2 = (r_dummy+skin)**2-r_dummy**2



!      
!---- Initialize random number generator
!
      iseed_z = iseed
      kptr = 1
      call inr250(mz,iseed,kptr)

!  ---- Read random number from file if exists
      
        iseed = iseed_z
        inquire(file='random_seed.dat',exist=es)
        if(es) then
            open(unit=222,file='random_seed.dat',status='old' )
            read(222,*) iseed
            close(222)
        end if
        print *,"  *  Initial random seed=",iseed
!      
! --- Initialize random number generator used for random force in DPD
!          (ziggurat)

#ifdef _OPENMP

allocate(par_jsrseed(0:numth-1))
call date_and_time(VALUES=values)
do j=0,numth-1
    par_jsrseed(j)=iseed+j+j*values(8)
end do 
call par_zigset( numth,par_jsrseed ,grainsize)
write(*,*) "  *  Parallel OpenMP version. N_threads: ", numth
#endif


      call zigset(iseed)

! Set zero the histogram vectors  

      histo_b(:) = 0.0 ; histo_d(:) = 0.0

! Random and dissipative forces for DPD

      force_r(:,:) = 0.
      force_d(:,:) = 0.


! NOTE: for DPD only friction(1) is used, the friction constant
! is assummed to be the same for all the particles

      sig = sqrt( 2.*temp*friction(1)/dt ) 

! NOTE: 1/sqrt(dt): is related to the random force in the algorhytm implemmentation 

! DPD constant weight functions. Here to speed up dpd_forces_ll

#if DPD_WEIGHT == 1
           w_d = 1.
           w_r = 1.
#endif   
  
!--------------------------------------------------
! Check consistency of using Poiseuille version of the program and moving the walls at the same time
#ifdef POISEUILLE
       if (va_spring_twall(1) /= 0._8 .or.  &
           va_spring_twall(2) /= 0._8 .or.  &
           va_spring_twall(3) /= 0._8 ) then
           print *, 'INCONSISTENCY: In Poiseuille flow version the walls should be at 0 velocity.'
           print *, 'Check mfa_input. Stopping here.'
           stop
       end if
#endif


! Parameters for molecules and beads

           part_init_d = n_mon*n_chain 
           part_init_e = n_mon*n_chain + n_mon_d*n_chain_d !NOTE: this is n_mon_tot if there is no particle 4
           part_init_star= part_init_e+ n_mon_e*n_chain_e

#if  SYSTEM == 1
           nm = n_mon_tot - part_init_d ! number of melt particles
           print '(/a,i6/)',"  * Melt particles= ",nm
#endif

! Params for observation

      if(mod(n_chain,2) == 0 ) then 
            n_chain_2 = n_chain/2
      else
        print '(/a/)', "WARNING: n_chain is not dividible by 2 ! "
        stop
      end if


  100 format(a)
  101 format(e14.6)
  103 format(3e14.6)
  110 format(i14)
  130 format(3i14)
  201 format(e14.6,a)
  210 format(i14,a)
  220 format(2i14,a)

      end subroutine init_params
