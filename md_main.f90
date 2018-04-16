program md_pb 
#include 'control_simulation.h'
    use commons 
    use util  

    implicit none
    integer  :: tot_time
    real (kind=8) :: t0,t1
    logical,parameter :: debug =.false.
    !!!DEBUG
    real*8, ALLOCATABLE :: f_old(:,:), f_new(:,:)
    integer :: i,j
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

!!! DEBUG

ALLOCATE(f_old(3,n_part))
ALLOCATE(f_new(3,n_part))

!!! DEBUG

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


        ! ----  Propagate coordinates 


        call verlet_positions()

        !deb           if (i_time == 1 ) then
        !deb               call write_conf(1,r0,10)
        !deb               stop
        !deb           end if


        call check_skin  !calculates if it is necessary to update the verlet list. If it is,
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


        ! NOTE: fluid-fluid calculates LJ always, DPD and LGV forces and calls ewald in
        ! real space
        call fluid_fluid() ! Calculates various forces, including thermostat and LJ


#if SYMMETRY == 0
#if WALL != 1 
        call fluid_wall(wall_flag) ! 1= wall atoms, 2= 9-3 potenti , 3 and 4 also valid
#endif

        call wall_wall(wall_flag)  ! 1= wall atoms, 2= 9-3 potential
#if WALL == 1
        call intra_wall
#endif
#endif
        call intra_molec


#ifdef ACTIVE_BRUSH
        call metronome(3) ! adds active brush forces
#endif

!!!DEBUG
f_old=force!(:,2+(n_chain-1)*n_mon)
!!!

#ifdef SPRING_ARRAY
        call spring_array(1) ! adds active brush forces
#endif

!!!DEBUG
f_new=force!(:,2+(n_chain-1)*n_mon)
!!!

#ifdef BENDING        
        call bending(1)  ! adds brush bending forces and bending energy
#endif


#ifdef BENDING_MELT        
        call bending_melt(1) ! adds melt bending forces and bending energy
#endif


#ifdef ORIENTATION        
        call orientation(1) ! adds brush orientation bending forces and bending energy
#endif


#if SYSTEM == 2 || SYSTEM == 3
        call ewald_k(1)  ! coulomb force calculation in K-space for Ewald sum
#if SYMMETRY == 0
        call dipolar_correction()
#endif
#endif

#ifdef POISEUILLE
        call constant_force() ! Poseuille flow generation
#endif


        ! -----  Update  velocities

        call verlet_velocities()
        !if(mod(i_time,100).eq.0) then
        !print *, a
        !endif
#ifdef DPD_VV                     

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


!TYPE vector
!    REAL :: x, y, z
!END TYPE vector

