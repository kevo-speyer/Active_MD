subroutine init_force_switch_on()
use commons 
implicit none
        real (kind=8) :: maxsig
! *** Setup 'force switch-on'

! r_2_min_init = effective mininum radius to be used for force computation at startup. This values disminishes
!                progressively with time in the variable r_2_min
! r_2_min_time = Maximum time until which the eff. radius will be used. 

!
      maxsig= maxval(sigma(:,:))

      r_2_min_init = (0.9*maxsig)**2 !ori= 0.9        !  Initial r_2 used for the  force "switch on" at first 
      print '(/a,f10.5/)',"  *  Effective min distance for 'force switch-on' = ",sqrt(r_2_min_init)

      if( dble(n_relax)*dt < 40.) then
          r_2_min_time = 40.               !  Real time after which the total force is completely switched on
      else    
          r_2_min_time = dble(n_relax)*dt   !  Real time after which the total force is switched on
      end if

      if ( int(r_2_min_time/dt) > n_relax+n_obser ) then
          print '(/a,f16.5,a,2x,i6//)',"  * WARN: The time to  &
              full force switch-on is higher than total simulation time."
      end if

      print '(a,f16.5,a,2x,i6//)',"  *  Force 'switching on' time= ",r_2_min_time," time steps= ",int(r_2_min_time/dt)

      ! *** Force "switch on" initialize ramp for effetive minimum distance

      !      r_2_min = r_2_min_init
  end subroutine init_force_switch_on
