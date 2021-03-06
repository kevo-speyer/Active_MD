program md_pb 
#include 'control_simulation.h'
    use commons 
    use util  

    implicit none
    integer  :: tot_time
    real (kind=8) :: t0,t1
    logical,parameter :: debug =.false.
    !!!DEBUG
!    real*8, ALLOCATABLE :: f_old(:,:), f_new(:,:)
!    integer :: i,j
    !!
    ! [12/2015] Paralelization with OpenMP
    ! [11/2013]: Droplet compatibility fixed
    ! [7/2013] : Bending Force for the brush added       
    ! March 2013: pass to Kevin
    ! [4/2010] Heavy changes in data structures and improvement of some routines 
    ! [5/2009]: Additions by Leonid Spirin: pinning, stars, LJ wall and shear
    ! protocols 
    ! [10/07/08]: Extended DPD cutoff, for better Temperature control.WARN: It changes
    ! a lot the friction 
    ! [20/11/07] Coulomb interactions added. It performs Ewald sum with dipole
    ! corrections to take into account the walls in Z. SYSTEM = 2 added.
    ! [18/11/07] Now with SYMMETRY = 1 the program has PBC in three directions for
    ! BULK simulations 
    ! [10/10/07] Particle 4 interaction added to the program. It can be also a  chain.
    ! [09/05/07] Fix an old error when writing confs. f_film discarded
    ! [2006] Incorporation of C-preprocessor directives
    ! [18/6/04] Pure velocity verlet  for integration and DPD thermostat implementation of
    !           Random and Dissipative forces
    ! [10/6/04] Shear: constant velocity for the heads of brushes. Force on heads re written.    
    ! [6/6/04]  Added a 9-3 flat wall potential (instead of explicit wall atoms) 
    ! [5/5/04]  Dynamic allocation of forces and position. New input file: system_input
    ! [14/4/04] claudio: using now ifort as the main compiler
    ! for development
    ! [then] Torsten Kreer : Brushes and melt with shear
    ! [originally from] Martin Muesser

    call get_walltime(t0)
    print '(/a/)'," *** MFA_PROG *** MD, DPD, and more. [Compilation Date: "//__DATE__//"]"


    ! ---- Write out current program compilation settings ----------     
    !
    call messages()

    ! ---- Setup some flags

    !      if(ad_flag.eq.0) then
    !       f_angle=0
    !      end if

    ! --------------- ENDS writing  out of compilation settings ---------------------------     

    ! ---- Inizialization  routines

    call init_system() ! physical system. Reads system_input
    call init_params()   ! md simulation parameters, box dimensions. Read mfa_input
    call init_config()   ! initial conf or read an old one
    call init_obser()    ! Measurable variables init.

!!!! DEBUG
!
!ALLOCATE(f_old(3,n_part))
!ALLOCATE(f_new(3,n_part))
!
!!!! DEBUG

#ifdef BENDING     
    call bending(0) ! writes brush bending constants  to log
#endif

#ifdef BENDING_MELT     
    call bending_melt(0) ! writes mel  bending constants  to log
#endif


#ifdef ORIENTATION     
    call orientation(0) ! writes brush bending orientation to log
#endif

#ifdef ACTIVE_BRUSH
    call metronome(0) ! add proper log info 
#endif

#ifdef SPRING_ARRAY
        call spring_array(0) ! adds active brush forces
#endif

    tot_time = n_relax + n_obser 

    !print *,"v_fluid",v_fluid_fluid ; stop

    do i_time = 1 , tot_time !  MAIN TIME LOOP 

        r_time = dble(s_time+i_time-1)*dt


        !DEBUG
        !  do i_part =1 , n_mon_tot
        !    print*,"force_long", force_long(:,i_part)
        !  end do
        !DEBUG
#ifdef RESPA
        call verlet_velocities_long(1)

        do i_time_short = 1, n_time_short ! inner (short) time loop for Multi-time-scale dynamics
       
#endif

        ! ----  Propagate coordinates 


        call verlet_positions() ! MAKE THIS RESPA COMPATIBLE!


        !!!DEBUG
        !do i_part = 1, n_mon_tot
        !print*,r0(:,i_part)
        !print*,v(:,i_part)
        !print*,a(:,i_part)+a_long(:,i_part)
        !end do
        !/DEBUG
        !deb           if (i_time == 1 ) then
        !deb               call write_conf(1,r0,10)
        !deb               stop
        !deb           end if


!DEBUG : HASTA ACA VA BIEN



#if BIN_TYPE == 0 || BIN_TYPE == 1
        call check_skin  !calculates if it is necessary to update the verlet list. If it is,
#endif

#if BIN_TYPE == 0

        if (f_skin.eq.1) call binning
#elif BIN_TYPE == 1 

        if (f_skin.eq.1) call my_binning
#endif


        ! ----- Forces to zero

        force(:,:) = 0.0 

        ! Note: After r_2_min_time, the force is completely switched on 
        ! This r_2_min is used for calculating the v_fluid_fluid and v_fluid_wall
        !
        ! Stuff for viscosity (it is spread in all routines that calculate forces and in
        ! viscosity.f90 )     ![VISC] 
        !
        !! WARN: We are assuming here that THIS routine is the first one in calculating
        !forces. That's why
        ! the pressure tensor is put to zero in each time step. IF this is not the first
        ! routine computing forces. 
        ! this should not be here.
#if SYMMETRY == 1
        press_tensor(:,:) = 0.
#endif

#ifdef RESPA
        !calc_short_force() calculates LJ between monomers-solvent & monomer-monomer and DPD or langevin
        call calc_short_force()
#else
        ! NOTE: fluid-fluid calculates LJ always, DPD and LGV forces and calls ewald in
        ! real space
        call fluid_fluid() ! Calculates various forces, including thermostat and LJ
#endif

#if SYMMETRY == 0
#   if WALL != 1
       call fluid_wall(wall_flag) ! 1= wall atoms, 2= 9-3 potenti , 3 and 4 also valid
#   endif

        call wall_wall(wall_flag)  ! 1= wall atoms, 2= 9-3 potential
#   if WALL == 1
#       ifdef RESPA
        print*,"ERROR: WALL = 1 not compatible with define RESPA."
        print*,"routine intra_wall should be compatibilized"
        stop
#       endif        
        call intra_wall
#   endif
#endif
       call intra_molec


#ifdef ACTIVE_BRUSH
call metronome(5) ! modify k_or, if Rend is in activation zone
#endif

!!!DEBUG
!f_old=force!(:,2+(n_chain-1)*n_mon)
!!!

#ifdef SPRING_ARRAY
       call spring_array(1) ! adds active brush forces
#endif

!!!DEBUG
!f_new=force!(:,2+(n_chain-1)*n_mon)
!!!

#ifdef BENDING       
       call bending(1)  ! adds brush bending forces and bending energy
#endif


#ifdef BENDING_MELT       
#   ifdef RESPA
    print*,"ERROR: BENDING_MELT NOT COMPATIBLE WITH DEFINED RESPA"
    stop
#   endif
        call bending_melt(1) ! adds melt bending forces and bending energy
#endif


#ifdef ORIENTATION      
       call orientation(1) ! adds brush orientation bending forces and bending energy
#endif


#if SYSTEM == 2 || SYSTEM == 3
#   ifdef RESPA
    print*,"ERROR: CHARGED SYSTEM NOT COMPATIBLE WITH DEFINED RESPA"
    stop
#   endif

        call ewald_k(1)  ! coulomb force calculation in K-space for Ewald sum
#   if SYMMETRY == 0
        call dipolar_correction()
#   endif
#endif

#ifdef POISEUILLE
       call constant_force() ! Poseuille flow generation
#endif
        !!DEBUG
        !do i_part = 1, n_mon_tot
        !print*,r0(:, i_part)
        !print*,v(:, i_part)
        !end do
        !!DEBUG

        ! -----  Update  velocities
        !MAKE THIS ROUTINE RESPA COMPATIBLE!
        call verlet_velocities()

#ifdef RESPA
        end do ! short time loop

        force_long(:,:) = 0.0 ! set solvent (long) forces to 0



!HASTA ACA VA BIEN 

        call calc_solv_solv_force() !
        
        !DEBUG
        !do i_part = 1, n_mon_tot
        !print*,force_long(:, i_part)
        !end do
        !DEBUG

        call verlet_velocities_long(2)
#endif        


#ifdef DPD_VV                     

#   if BIN_TYPE == 2
print "ERROR: DPD_VV not compatible with BIN_TYPE 2 yet!"
print*, "Exiting program"
exit
#   endif
        !Note: this recalculates Fd with the new velocities and updates F for the begining og the next cycle 
        !       with this new value.  

        call new_dpd_fd()  
#endif

        !----  Observe system after equilibration

        if(i_time.gt.n_relax) call observation 
        
        !----  Make safety copies  to recover from crashes and write out of configurations

        if(mod(i_time,n_safe).eq.0) then
            call store_config(2)  ! writes out conf_xmol and conf_new

        !!<--- DEBUG ----
        !open(143,file="force_xmol",status="unknown",position="append")
        !write(143,*) n_part
        !write(143,*) 
        !do j=1,n_part
        !    write(143,'(A)',advance='no') 'Type'
        !    do i=1,3
        !        write(143,203,advance='no') f_old(i,j)-f_new(i,j)!'(ES12.4)' 
        !    end do
        !    write(143,*) 
        !end do
        !close(143)
        !!---- DEBUG --->


#if STORE == 0
            call store_config(3)  ! Writes out film_xmol and vel.dat
#elif STORE == 1
            call store_config(4)  ! Writes out film_xmol and vel.dat UNFOLDED
#endif
            call synchro_obser(2)

        end if

    end do   ! --------------  ENDS TIME LOOP ----------------

    call obser_out()

    close(20) ! closing mfa_output
    !
    call get_walltime(t1)
    print *,'    WALL TIME (s)= ',t1-t0

    203 format(3f13.4)
end program md_pb

