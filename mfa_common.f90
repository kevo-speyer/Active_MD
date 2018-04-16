module commons
      implicit none
      save ! is not strictly necesary because commons will be used also in the main program 

#include 'control_simulation.h'
      
!---  conventions about variables:
!---------------------------------


!---  no variables with just one letter 
!     (hard to find in editor and hard to replace as well)

!---  variables starting with: f_
!               flaggs

!---  variables starting with: n_
!               dimensions of arrays & lenghts of loops

!---  variables starting with: i_  j_  k_ l_
!               loop indices (not allowed in! common blocks)
!               not to be passed on to subroutines (might be changed)

!---  variables ending with: _1 _2 _3 _4
!               1., 2., ... power or cumulant of a variable

!---  pass variables to subroutines by call <subr>(x1,x2,..)
!     only if subroutine <subr> doesn't read the! common block file



!---  conventions about files:
!-----------------------------


!---  input files:  name <file>.old
!     the first number in an input file is either zero or one
!     zero: no information is contained and default values are taken
!     one:  complete information is available which will be read in

! Ordering of variables, to ease cache alignement

! ---- INTEGER PARAMETERS

      integer, parameter  ::  n_dim=3

#ifndef STARS
#   if WALL==1
#        if SYSTEM ==3 || SYSTEM==0 || SYSTEM==1 
#           ifdef  PARTICLE_4
              integer, parameter  :: n_type=5
#           else               
              integer, parameter  :: n_type=4
#           endif
#        endif
#   elif WALL==2 || WALL==3 /* if wall is implicit  */
              integer, parameter  :: n_type=4
#   endif
#else /* if there are stars */
         integer, parameter  :: n_type=7
#endif
#if SYSTEM == 1 
!              integer, parameter  :: n_type=4
#endif              

#if PINNED==1 || PINNED==2
      real(kind= 8) ::  dist_pinned, r_start(2, n_dim)
#endif

      integer, parameter  :: n_order=2
      integer, parameter  :: n_layer =1

      integer, parameter :: n_diff_time=9    

      integer, parameter :: n_bin_wa=16,            & !ori: 16
                            n_bin_fl= 300    !ori: 256       !increíble: con 310 ya no compila                     

      integer, parameter :: n_neigh_fl=256,        &  !n_neigh_fl=256       !ori
                            n_neigh_wa=32,         &  !n_neigh_wa=32,       !ori
                            n_neigh_ww=32             !n_neigh_ww=32        !ori

! Histogram variables for brushes (_b) and droplet/melt (_d)

      integer, parameter :: hist_dim = 200

      integer,parameter ::  n_layers = 200

      integer, parameter         :: wall_flag = WALL   
! Paramas for optimizition

!      integer, parameter :: csize = BLOCK_SIZE

!
! DOUBLE PRECISION VARIABLES 
!
 

      real (kind=8) ,allocatable ,target :: r0(:,:),spabs(:,:)
#ifdef STARS
      
      real (kind=8) ,allocatable  :: r0_star(:,:)
#endif
      real (kind=8) ,allocatable  :: r0_unfold(:,:)

      real(kind=8), allocatable ::  r_wall_equi(:,:)


      real (kind=8) , allocatable :: r_head_old(:,:)                          

       real (kind=8) ,allocatable :: v(:,:),a(:,:),v_half(:,:)
#ifdef PROFILES
       real (kind=8) ,allocatable :: v_prof(:,:,:),v_prof_2(:,:,:),histo_out(:,:)
#endif

! force on head brushes to calculate shear forces
      real (kind=8), allocatable :: f_on_heads(:,:)


      real (kind=8)     ::   r_time
      real (kind=8)     :: drop_cell_x,drop_cell_y 

!---  definitions for intra-wall spacing
      real (kind=8)     ::  x_space,y_space,x_max,y_max,surface
#if SYSTEM == 1
      real(kind=8) :: x_max_d,y_max_d  ! DROP

#endif /*droplets*/

!---  definitions for inter-wall spacing
      real (kind=8) :: z_space_wall
      real (kind=8) :: z_head ! Z position of the brush heads

      real (kind=8) ::  v_wall_wall_1,t_wall_1,v_wall_wall_2,t_wall_2,e_wall_wall_1,e_wall_wall_2
!---  variables needed for per.bound.cond. and computation of potential
      real (kind=8) :: boundary(3),inv_boundary(3),half_boundary(3)
#if SYSTEM==1 /*DROPLETS*/
      real (kind=8) :: boundary_d(3)
      real (kind=8) :: z_max_d,z_min_d,x_box_min,x_box_max,vol_drop,dens_drop
#endif
!---  fluid wall observables
      real (kind=8) :: v_fluid_wall,v_fluid_wall_1,v_fluid_wall_2,inv_N,inv_V,inv_n_mon,inv_n_chain
!---  fluid fluid observables
      real(kind=8)  ::  v_fluid_fluid,v_fluid_fluid_1,v_fluid_fluid_2,t_fluid, &
                        t_fluid_1,t_fluid_2,e_fluid_fluid_1,e_fluid_fluid_2
!---  global observables
      real (kind=8) ::  v_total,v_total_1,v_total_2,t_total,t_total_1, t_total_2,e_total,e_total_1,e_total_2
      real (kind=8) ::  v_tot_intra_1,v_tot_intra_2,v_tot_intra_d_1,v_tot_intra_d_2

!---  variables for top wall
      real (kind=8) :: r0_twall(n_dim)
!---  variables for bottom wall
      real (kind=8) ::  r0_bwall(n_dim)

      real (kind=8) :: r0_spring_twall(n_dim)

!! obs ?      real (kind=8) :: diff_lev1(n_dim,n_diff_time),diff_lev2(n_dim,n_diff_time), &
!! obs ?                       diff_lev3(n_dim,n_diff_time),diff_lev4(n_dim,n_diff_time), &
!! obs ?                       diff_lev5(n_dim,n_diff_time),diff_lev6(n_dim,n_diff_time)
! skin of the binnig boxes
      real (kind=8) :: binning_box_skin    

!   Variables for my_binning 
    integer :: n_cell 
    integer,allocatable :: bin(:),binpnt(:)
!obs?      real (kind=8) :: currdens(0:11,0:2),currdensup(0:11,0:2),currdensdown(0:11,0:2),dens(1:11),densup(1:11),densdown(1:11)

! parameters for the wall-fluid interaction and potential.
! Used when no explicit wall particles defined in the model
#if WALL==2 || WALL == 3
      real (kind=8) :: sigma_w,a_w,sigma_w4,a_w4, sigma_wall(n_type), a_wall(n_type)
#endif

      real (kind=8) :: histo_b(hist_dim),histo_d(hist_dim)
! Rcut variable
      real (kind=8) :: r_cut_max,vec_dummy(3),z_skin,r_cut_dpd,r_cut_max_2,r_cut_dpd_2
!$OMP THREADPRIVATE(vec_dummy)      
! Radius of gyration parameters            !(r_g² or other quantity, part_type,coor)
      real (kind=8) :: r_cm(3),r_g2(2,3),r_g2_mean(3,3),r_par_per2(3),v_dummy(3),re_2_mean(3,3)

      real (kind=8) :: mean_f_on_heads(3)

! DPD implementation with Velocity Verlet integration 

! Metronome routine variables
REAL(8) :: pot
REAL (KIND=8) :: k_act
# ifdef ACTIVE_BRUSH
REAL (KIND=8), ALLOCATABLE ::  cos_mem(:)
INTEGER, ALLOCATABLE ::  kick_mem(:)
REAL (KIND=8), ALLOCATABLE ::  k0(:)
# endif

REAL (KIND=8) :: r_rel(3), cos_th ! For angle monitoring in fort.56

! Spring array routine variables
REAL (KIND=8) :: k_spr_x, k_spr_y, v_array ! spring elastic constants
# ifdef SPRING_ARRAY
INTEGER, ALLOCATABLE :: Mindex(:,:) ! Index matrix of 2 bead of i-chain for r0
# endif

! Bending parameters and forces
# ifdef BENDING        
real (kind=8) :: k_bend, alpha_eq, v_bend,v_bend_melt  !elastic constant, equilibrium angle
real (kind=8) ,allocatable :: force_bend(:,:)
real (kind=8) :: beta1, beta2 ! angles for active brush model.
# endif
# ifdef ORIENTATION        
real (kind=8) :: k_or, alpha_or,v_or  !orientation elastic constant, orientation equilibrium angle, orientation potential energy
real (kind=8) ,allocatable :: force_or(:,:)
# endif


! Vectors to store the Random and Dissipative force

      real (kind=8), allocatable :: force_d(:,:), force_r(:,:)
!obs      real(kind=8), allocatable  :: random_2(:),random_3(:)

! DPD weight funtions variables 

       real (kind=8) :: w_d, w_r
!$OMP THREADPRIVATE(w_d,w_r)      

! Coeficients for the new integrator 

       real(kind=8) :: sig 

       real(kind=8) :: tot_P(3),tot_P_2(3)
       real (kind=8) :: D_coeff,D_coeff_cm

#ifdef PROFILES
! Put here all that related to mide for calculating profiles 
       real (kind=8) :: zz_min
#endif

! POISEUILLE
! Variables for the constant external force simulation [MFA_POISS ]
       real(kind=8) :: const_force(3) ! It is physically a gravity. Equal to the force if all the masses are identical
       real(kind=8) , allocatable :: ext_force(:,:)


! Number of layers in the different profiles (velocity, etc)
! ------  REALS


!---  veloc., accel. and higher derivatives

!      real (kind=8), allocatable :: rx(:,:,:)
      real (kind=8), allocatable :: force(:,:)
      real (kind=8), allocatable :: mass(:),inv_mass(:)

! This is OK with real (not double precision!) 
      real, allocatable :: random(:)

!      real, allocatable :: r_level1(:,:,:),                    &
!                     r_level2(:,:,:),                    &
!                     r_level3(:,:,:),                    &
!                     r_level4(:,:,:),                    &
!                     r_level5(:,:,:),                    &
!                     r_level6(:,:,:)

      real(kind=8) ,allocatable:: r0_old(:,:)
!---  variables for tensor of gyration
!obs      real(kind=8), allocatable ::  chcm(:,:), gyrten(:,:,:), gyrt(:,:,:)   

! Diffusion coeficient variables 

#ifdef DIFF
      real(kind=8),  allocatable ::  r0_ref(:,:)  ! ,r0_cm_m(:,:)
#endif      

      real (kind=8) :: mass_type(n_type),    &
      epsil(n_type,n_type),sigma(n_type,n_type),range_2(n_type,n_type),e_shift(n_type,n_type), &
              friction(n_type),rf_width(n_type),sigma_2(n_type,n_type)
#       if SOLVENT == 2 || SOLVENT == 3 
! variables for chemical incompatibility and non-additive potentials 

        real (kind=8) :: delta_sig(3)
        
#       endif
      real (kind=8) :: inv_range_2(n_type,n_type)

!---  predictor corrector coefficients
      real predict_coef(0:n_order-1,0:n_order),correct_coef(0:n_order)

      real(kind=8)   ::  dt,dt_2,inv_dt,inv_dt_2
      real(kind=8)   ::  r_dummy,r_dummy1,r_dummy2
      real(kind=8)   ::  k_spring_wall
!---  wall observation variables
      real(kind=8)   ::  v_wall_wall,t_wall

!---  intramolecular properties
      real (kind=8)   :: r_chain,k_chain,r_chain_2,v_intra_molec,v_intra_molec_d,v_intra_molec_e,inv_r_chain_2
!---  ramp for interaction
      real(kind=8)   :: r_2_min_init, r_2_min_time, r_2_min !(n_type, n_type)
!---  temperature
      real(kind=8)   :: temp !obsoleted ,temp_time,temp_init,temp_final

       real (kind=8)  :: delta_r(n_dim),force_loc(n_dim),r_2,r_6,r_12,pot_loc
!$OMP THREADPRIVATE(r_2,delta_r)      


      real(kind=8)            :: rx_twall(n_dim,n_order),mass_twall
      real(kind=8)            :: va_spring_twall(n_dim),k_spring_twall(n_dim), force_twall(n_dim),r0_pinned_twall(n_dim)

      real(kind=8)            :: r0_last_pinned(n_dim),fo_last_pinned(n_dim) ,fo_last_unpinned(n_dim)
      real(kind=8)            :: r0_twall_1(n_dim),f0_spring_1(n_dim)

      real(kind=8)            :: r_bin_x,r_bin_y,r_bin_z   ! binning box dimensions
      real(kind=8)            :: inv_r_bin_x,inv_r_bin_y,inv_r_bin_z   ! HPC

      real(kind=8)            :: skin,skin_2

      real(kind=8)            :: z_1,z_2

      real(kind=8)            :: binbox,binstep,i_histodown(0:110),i_histoup(0:110)
!---  variables for force on top wall
      real(kind=8)            :: ftw(1:n_dim)
! OBSOLETE !---  variables for force from down brush on upper brush
! OBSOLETE       real            :: fordbub(1:n_dim)

      real(kind=8)            :: binbox_2,binstep_2 
      real(kind=8)            :: crdna(1:3),crdnb(1:3),rend(1:3),rendtot,angle, rendl2,zcomp,xcomp,zcoord
      real(kind=8)            :: crdna2(1:3),crdnb2(1:3),rend2(1:3),rendtot2,angle2, rendl22,zcomp2,xcomp2,zcoord2

! Some variables for measuring quantities

      real(kind=8)            :: v4_mean(3),v4_mean2(3)
      integer                 :: c4
!
! ---- REAL PARAMETERS
! 

      real(kind=8), parameter :: pi=3.1415926536,twopi=2*pi,sqrt2=1.41421356237,sqrt3=1.73205080756
!
! ------ INTEGERS 
!

!---  variables to monitor mimimum image convention
     integer,allocatable ::  mic_count(:,:)

!---  definitions for particles
      integer i_type,j_type 
!$OMP THREADPRIVATE(i_type,j_type)      
      integer,allocatable  :: a_type(:)
!--- Binning params
      integer , allocatable ::  bin_twall(:,:,:,:)
      integer , allocatable ::  bin_bwall(:,:,:,:)
      integer , allocatable ::  bin_fluid(:,:,:,:)       

      integer, allocatable  ::   ff_list(:,:) ! ,                         &
!                                fw_list(:,:),                         &
!                                ww_list(:,:)


      integer, allocatable  :: mic_old(:,:)

      integer :: i_dummy, j_dummy

!---  s_time:  starting integer time; n_relax: number of relax. steps
!     n_obser: number of steps in which observation takes place
!     r_time = i_time * dt 
      integer :: i_time,s_time,n_relax,n_obser !,counttime

! System characteristics
! ---------------------
!---  enumeration of monomers (fluid), chains, wall atoms, ...
!     n_mon: number of monomers on one chain
!     n_chain: number of chains
!     n_wall: number of wall atoms
!     n_cell_w_x,n_cell_w_y : number of cells for wall in x,y direction

      integer :: i_mon,i_chain,i_mon_tot ! ,n_mon_tot,n_mon,n_chain  !ori
      integer :: i_cell_w_x,i_cell_w_y    !,n_cell_w_x,n_cell_w_y !ori
      integer :: i_wall,i_part,j_part,q_part !,n_part,n_layer !,n_wall !ori
!$OMP THREADPRIVATE(j_part,q_part)     
      integer :: n_safe
!obs      integer :: f_wall_fix,f_sp_wa,f_th_fl,f_th_wa,f_up_li,f_film, ad_flag,f_profile,f_force,f_curr,f_logscale,f_angle 
      !integer ::  !ad_flag !,f_wall_fix !,f_minimize
!OBSOLETE: solv_flag
! OBSOLETE      integer :: f_wiw,f_intcount,g_gyrten 

      integer ::            n_cell_w_x,n_cell_w_y, n_wall  

! Definition of limit variables to be read or calculated in init_system
      integer :: n_mon,n_chain,n_mon_tot,n_part                    

! Definition of  droplet variables: (or the melt according the shape)
      integer :: n_mon_d,n_chain_d,part_init_d,n_chain_e,n_mon_e,part_init_e,n_tot_e, part_init_star
      integer :: n_chain_2
#ifdef SHEARED 
      integer ::  n_mon_arm, n_arms, n_stars, i_star, i_arm, i_mon_arm, d_part
      integer :: i_mon_d, i_chain_d
      integer :: turn_time(3), p_time

#endif
!Parameters for shear Protocols
      real (kind=8) appl_vel(n_dim),appl_vel_bw(n_dim),appl_vel_tw(n_dim)
      real (kind=8) appl_acc(n_dim),appl_acc_bw(n_dim),appl_acc_tw(n_dim)
      real (kind=8) shifttw(n_dim),shiftbw(n_dim), fbw(1:n_dim)

#ifdef PARTICLE_4
! Definition of particle type 4 if used
!           Chain length, number of chains and type of interaction of particle 4
!  1= good solvent, i.e repulsion with cut off 1.2 sigma with all the other
!  types. 
!
!  - Defined in init_params.f90 

      integer :: inter4
#endif
!      integer           ::  n_loop
      integer           ::  i_dim,j_dim
      integer           :: i_order,j_order 
!obsolete      integer           :: f_cut_off
      integer           :: iseed,mz(250),kptr
      integer           :: count_obs
      integer           :: pbc_twall(n_dim)
      integer           :: pbc_bwall(n_dim)

!---  variables for external side conditions, including external spring
      integer           :: f_twall(n_dim)

      integer           :: s_force_grad(n_dim)

!---  time period over which trajectory is averaged
      integer           :: time_ave_count,time_ave_count_2,n_time_ave,n_time_ave2

!---  variables to monitor diffusion ![ori not used/understood ]

!obs      integer           :: i_diff_time
!obs      integer           :: i_obs_l1,i_obs_l2,i_obs_l3,i_obs_l4,i_obs_l5,i_obs_l6

!---  variables for binning
      integer           :: i_bin_x,j_bin_x,k_bin_x
      integer           :: i_bin_y,j_bin_y,k_bin_y 
      integer           :: i_bin_z,j_bin_z,k_bin_z  
      integer           ::  n_bin_x, n_bin_y,n_bin_z                         

      integer  :: i_bin_fl,j_bin_fl   ! i_bin_wa,n_bin_fl,n_bin_wa         !ori

      integer            :: i1_bin_part,i2_bin_part,n1_bin_part,n2_bin_part
      integer            :: delta_bin_x,delta_bin_y,delta_bin_z,n_fcell_max,n_wcell_max
!---  variables for neighbor list
      integer            :: i_neigh, n_neigh_fl_max,n_neigh_fw_max,n_neigh_ww_max
      !integer            :: n_list_up

       integer           :: f_skin

!---  variables for monomer profile
      !integer  :: n_step !,i_step

!OBSOLETE !---  variables for interaction counting
!OBSOLETE       integer            :: intcount      
!---  variables for current density and related density profile
!obs      integer            :: n_step_2

      integer            :: bingyr,countgyr(1:10) ! radius of gyr variables 

!---  variables for angle of tilting and end-end-vector
!obs      integer            :: cntr,cntr2
!obs      integer            :: tcount
    
                           ! according to f_explicit_wall
#ifdef PROFILES
      integer            :: dim_prof_dens
#endif

!--------------------------------------------------------------------------
! NEW SWITCHS. Change substantially the program scheme !
! Parameter to  control explicit wall atom simulation or 9-3 pot or explicit wall
! true =  explicit wall  atoms 
! false= 9-3 potential, flat wall

!!! obsolete #if WALL == 2 || WALL == 3 || WALL == 4
!!! obsolete         logical,parameter :: f_explicit_wall=.false.    ! *** important switch ***   
!!! obsolete #endif
!!! obsolete #if WALL == 1
!!! obsolete         logical,parameter :: f_explicit_wall=.true.    ! *** important switch ***   
!!! obsolete #endif

! false=  original integration scheme and thermostat of the program 
!         (predictor corrector max_order=5 n_order=2 ~velocity Verlet) 
! true=   DPD thermostat and different integration scheme ~Velocity Verlet but
!         friction and random terms are updated at the positions stage
!!! #if THERMOSTAT == 0 
!!!         logical,parameter :: f_dpd=.true.              ! *** important switch ***   
!!! #endif        
!!! #if THERMOSTAT == 1 
!!!         logical,parameter :: f_dpd=.false.              ! *** important switch ***   
!!! #endif        
!! obs #if THERMOSTAT == 1 && defined LGVX_0
!! obs ! WARN: puts 0 friction force in the X direction !!!!!!!!!!!!!        
!! obs ! LANGEVIN is applied only in y and z directions  
!! obs         logical,parameter :: f_g0=.true.              ! *** important switch ***   
!! obs #else
!! obs         logical,parameter :: f_g0=.false.              ! *** important switch ***   
!! obs #endif        

! ---- CHARACTERS 

      character (len=64) ::  c_dummy                 
!---  special variables
      character (len=20) ::  name

!---  flaggs: 0 (off)   1 (on)
!     f_wall_fix: wall atoms constrained (1) or vibrating (0)
!     f_sp_wa: coupling from spring to wall
!     f_th_fl: thermostat fluid
!     f_th_wa: thermostat wall
!     f_up_li: update lists
!     f_film:  generate film_xmol (1)
!     f_wiw: distinguish upper and lower brush
!     f_profile: calculate monomer profile (1)
!     f_force: write out forces acting on top wall (1)
!     f_curr: calculate current dansity (1)
!     f_intcount: count number of interactions between brushes   OBSOLETE
!     f_gyrten: measure tensor of gyration
!     f_logscale: write out on logscale 
!     ad_flag: polymers adsorbed (0) or endgrafted (1)
!OBSOLETE     solv_flag: good (0) or bad (1) solvent without explizit
!OBSOLETE                fluid particles
!     f_angle: calculate mean angle of tilting (1)
!      common /cflaggs/ f_wall_fix,f_sp_wa,f_th_fl, &
!     &                 f_th_wa,f_up_li,f_film,ad_flag,solv_flag, &
!     &f_profile,f_wiw,f_force,f_curr,f_intcount,f_gyrten,f_logscale, &
!     &f_angle
 

#if SYSTEM == 2  || SYSTEM == 3
! Coulomb variables definitions
        integer ::  n_kmax
        real(kind=8), allocatable  :: q(:)
        real(kind=8) :: v_coul,v_coul_self,alpha_c,sqrt_alpha_c,charge(n_type) ,inv_sqrt_pi
        real(kind=8) :: sqrt_alpha,lz_half,z_slab
#endif

#if SYSTEM == 1 
! ---- Droplet definitions
! Vapour phase variables. This is droplet + low dens. vapour coexistence

      real (kind=8) :: vap_dens,vap_vol,vap_vol_box_1
      integer :: n_mon_tot_v,n_mon_vap,n_chain_vap


! Variables for measuring DROPLET characteristics and profiles

      real(kind=8) :: vcm_d(3), vcm_d2(3)

! Taken from mide.f90 from mide_drop_force

      integer :: n_free_chain,n_box(3),n_box_xz(3) ! may be not used. Placed here to
      
      real (kind=8) :: binning_box_length,guessed_drop_cm(3),r_box_xz(3),Lx_o_2(3)
      real (kind=8),allocatable :: dens_xz(:,:,:),dens_xz_b(:,:,:),dens_xz_2(:,:,:),dens_xz_b_2(:,:,:) 

!NOT YET #ifdef FORCE_PROFILE
!NOT YET real (kind=8),allocatable :: force_prof(:,:)
!NOT YET integer, allocatable :: histo_f(:)
!NOT YET real(kind=8) :: fix_shift
!NOT YET #endif
      real (kind=8),allocatable :: vel_xz(:,:,:),vel_xz_b(:,:,:), vel_xz2(:,:,:),vel_xz_b2(:,:,:)

! used in rebuild_drop2 or rebuild_drop
      real (kind=8) :: drop_cm(3),x_droplet_shift,fix_shift
      integer :: n_cm,nm

#endif /*SYSTEM=1: droplets*/

      real (kind=8), allocatable :: r_cms(:,:)
      
#if SYMMETRY == 1 
    
! Difussion coeficient variables 
!not yet logical :: ref_file, unfold_file  ! for difussion we need global acces of this 
!not yet real,  allocatable ::  r0_ref(:,:),r0_cm_m(:,:),r0_unfold(:,:)
!not yet real (kind=8) :: D_coeff,D_coeff_cm

! Viscosity calculatoin variables
! Pressure tensor = - shear tensor
! Must be global because the components ared added in each routine that calculates forces. 

real (kind=8) :: press_tensor(3,3),press_tensor_mean(3,3) 


#endif

!!OMP VARIABLES 
#ifdef _OPENMP
integer, dimension(:), allocatable :: par_jsrseed
integer :: ith, numth, grainsize=64
!$OMP THREADPRIVATE(ith)     
integer,dimension(8) :: values
#endif

      end module commons
