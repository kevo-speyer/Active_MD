
subroutine messages 
! * This routine writes to standard output the compilation settings of the program
! based on the variables defined in control_simulation.h
! * Performs consitency checks at compilation time.  
      use commons 
      implicit none
      logical,parameter :: debug =.false.

#include 'control_simulation.h'

! Symmetry of the problem: PBC in 3D or 2D

#if SYMMETRY == 0 
        print '(/a/)', "  *  SIMULATING: a system with PBC in X and Y coordinates (Walls in Z)" 
#elif SYMMETRY == 1
        print '(/a/)', "  *  SIMULATING: a system with PBC in X, Y and Z coordinates (BULK)" 
#endif

! --- Artifical settings ---- 

#ifdef BIDIMENSIONAL
        print '(/a/)', "  *  SIMULATING: Velocities in Y coordinate are permanently set to ZERO, and won't evolve"
#endif

! --- System settings ---- 

#if SYSTEM == 0
#   if SYMMETRY == 0 
        print '(/a)', "  *  SIMULATING: grafted chains, non grafted chains and 3 or 4 type of particles"
        print '(a/)', "  *  Channel-like geometry of the simulation box"
#   else
# ifdef ACTIVE_BRUSH
        print '(/a)', "  *  SIMULATING: active, polymers are driven by external forces"
# endif
# ifdef SPRING_ARRAY
        print '(/a)', "  *  SIMULATING: springs, brush polymers are joined by springs. Only works for ordered brush "
#endif
#   endif

#elif SYSTEM == 1 
        print '(/a)', "  *  SIMULATING: grafted chains, non grafted chains and 3 or  type of particles"
        print '(a/)', "  *  Droplet geometry: No brush in the top wall."
        print '(a/)', "  *  Droplet geometry: Others and more measuring routines"
#elif SYSTEM == 2
        print '(/a)', "  *  SIMULATING: free chains, ions and counterions with Coulomb interactions"
         print '(a)', "  *  No brushes."

#       if SYMMETRY == 1 
        print '(a/)', "  *  Performing Ewald sum in 3D" 
#       elif SYMMETRY == 0 
        print '(a)',  "  *  Performing Ewald sum in 3D"
        print '(a/)', "  *  Using dipolar correction in Z direction"
#       endif
#elif SYSTEM == 3
        print '(/a)', "  *  SIMULATING: brushes,   free chains, ions and counterions with Coulomb interactions"
        print '(a/)', "  *  Channel-like geometry of the simulation box"
#       if SYMMETRY == 1 
        print '(a/)', "  * bulk + fixed brushes is not possible ! CHANGE control_simulation.h !!! "  
        stop 
#       elif SYMMETRY == 0 
        print '(a)',  "  *  Performing Ewald sum in 3D"
        print '(a/)', "  *  Using dipolar correction in Z direction"
#       endif

#endif

! Messages related to the wall interactions
#   if SYMMETRY == 0
#       if WALL == 1 
               print '(/a/)', "  *  WALLS: Explicit LJ particles "
#       endif
#       if WALL == 2 
               print '(/a/)', "  *  WALLS: Implicit 9-3 potential. It could be purely repulsive if A_w<0"
#       endif
#       if WALL == 3 
               print '(/a/)', "  *  WALLS: Bottom: Implicit 9-3 potential. TOP: hard wall (bounce in z) "
#       endif
#       if WALL == 4 
               print '(/a)', "  *  WALLS: Both, TOP and BOTTOM walls are hard walls."
               print '(a/)', "  *  A_w and s_w are ignored. "
#       endif
#   endif   


#ifdef RELAX
        print '(/a)', "  *  WARNING:  the program is compiled in RELAX mode. NO MD !!! "
        print '(a)', "  *  WARNING:  the program is compiled in RELAX mode. NO MD !!! "
        print '(a/)', "  *  WARNING:  the program is compiled in RELAX mode. NO MD !!! "
#endif
#if SOLVENT == 0      
      print '(/a/)',"  * Simulating poor solvent conditions "

#elif SOLVENT == 1
      print '(/a/)',"  * Simulating good solvent conditions "
#elif SOLVENT == 2
      print '(/a/)',"  * Simulating good solvent conditions with non-additive potentials. Different sig_23  "
#elif SOLVENT == 3
      print '(/a)',"  * Simulating poor solvent conditions with non-additive potentials. Different sig_23  "
!      print '(a/)',"  * Interaction 2-3 is purely repulsive "
#endif
#ifdef PROFILES 
print '(/a/)',"  * Doing simulation with profiles calculation (slower) "
#endif

#ifdef FENE_l0
    print '(a)',"  * Enforcing minimum distance in FENE interactions"
#endif

#if THERMOSTAT == 0
      print '(/a)',"  * Using DPD thermostat. "

#   ifdef DPD_VV
      print '(a)',"  * Using DPD_VV: Fd is recalculated after the standard Velocity Verlet loop"
#   endif      
#   ifdef DPD_CUT_OFF
      print '(/a,f10.5/ )',"  * Using DPD cutoff independently of potentials: dpd_cut_off= ",DPD_CUT_OFF 
#   else      
      print '(/a/ )',"  * Using DPD cutoff coincides with potentials cut-offs" 
#   endif      

! ---- Different weight functions      

#   if DPD_WEIGHT == 0 
      print '(a)',"  * Using usual DPD weight function: Wd=(1-r/rc)^2 ; Wr^2=Wd "
#   endif      
#   if DPD_WEIGHT == 1 
      print '(a)',"  * Using CONSTANT DPD weight function: Wd=1 ; Wr^2=Wd "
#   endif      
#   if DPD_WEIGHT ==2 
      print '(a)',"  * Using SQUARE_ROOT  DPD weight function: Wd=sqrt(1-r/rc) ; Wr^2=Wd "
#   endif      
#   if DPD_WEIGHT ==3 
      print '(a)',"  * Using POWER 1/4  DPD weight function: Wd=(1-r/rc)^(1/4) ; Wr^2=Wd "
#   endif      
#   if DPD_WEIGHT ==4 
      print '(a)',"  * Using POWER 1/4*GAUSSIAN DPD weight function: Wd=(1-r/rc)^(1/4)*2.6*exp(-(1-r/rc)^2) ; Wr^2=Wd "
#   endif      
#   ifndef DPD_WEIGHT
      print '(a)',"  * Current version of the programs NEEDS a value for DPD_WEIGHT in dpd mode."
      stop
#   endif      
!!!! #ifdef DPD_EMBEDDED
!!!!         print '(a/)',"  * DPD forces calculations done INSIDE the fluid-fluid routine"
!!!! #else
!!!!         print '(a/)',"  * DPD forces calculation is done in dpd_forces_ll.f90"
!!!! #endif

#elif THERMOSTAT == 1
     print '(/a/)',"  * Using LANGEVIN thermostat. "
#endif

#if SYSTEM == 0 
#   ifdef POISEUILLE
           print '(/a/)',"  * POISEUILLE version of the program. External force used and expected in system_input"
#   else 
#       if SYMMETRY != 1           
           print '(/a/)',"  * Shear version of the program (Couette). Walls can be moved at V=const."
#       endif           
#   endif
#endif

#ifdef NO_WARNS
        print '(/a/)',"  * No warnings will be sent to standard output" 
#else
        print '(/a/)',"  * Some warnings will be sent to standard output" 
#endif
!! #ifndef FLUID_ROUTINE
!!            print '(/a/)' , "FLUID_ROUTINE not defined !!! Stopping here ! " ; stop
!! #endif
!! #if FLUID_ROUTINE == 0 
!!            print '(/a/)' , "  *  Using stable fluid_fluid.f90 version" 
!! #elif  FLUID_ROUTINE == 1
!!            print '(/a/)' , "  *  Using TESTING fluid_fluid_test.f90 version" 
!! #endif

#if STORE == 0 
           print '(/a/)' , "  *  Writing out folded coordinates"
#elif STORE == 1
           print '(/a/)' , "  *  Writing out UNFOLDED coordinates. This reads/writes conf_unfold"
#     ifdef DIFF
           print '(/a/)' , "  *  Performing diffusion calculation of particle 4: Each run starts a new reference "
#     endif
#endif

!! if(f_g0 ) then 
!!     print '(/a/)', "  * FRICTION force in X direction (shear) is set to 0 "  
!!     print '(/a/)', "  ---- IMPORTANT WHEN USING LANGEVIN THERMOSTAT !!!   ---- "
!! end if

! Consistency checks of the flags defined in control_simulation.h for  a given
! physics. 

#ifdef STARS
#    ifdef FLUKT
      write(*,*) "can not be Fluktuations with stars, not in program yet,"
      write (*,*) "to add it - just the same like in balls, it is just writing out fort.555"
     stop
#    endif
#    if SYSTEM == 1 || SYSTEM == 2
     write(*,*) "SYSTEM 2 or 3 with stars are not possible yet"
     stop
#    endif 
#    if WALL==2
     write(*,*) "Stars are only with explicit walls, for implicit wall modify types, and check reading of parameters from mfa_input." 
     stop
#    endif 
#endif

#ifdef FORCE_SWITCH_ON
        print '(/a/)',"  *  Performing progressive force switch on towards thermalization  " 
#endif

#if THERMOSTAT == 1 
!!!! #ifdef DPD_EMBEDDED
!!!!         print '(/a/)', "DEFINITION ERROR: with LGV thermostat  DPD_EMBEDDED must be undefined"
!!!!         stop
!!!! #endif
#endif

#if WALL==1
#    if SYSTEM==2 || SYSTEM==1

            print *,"  **** Incorrect algorithm combination: ***  "
            print *," explicit walls atoms with  system 2 or system 1 are not implemented ! "
            print *, "Stopping here. "
           stop

#    endif
#endif







end subroutine messages
