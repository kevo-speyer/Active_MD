      subroutine store_config(mode)
      use commons ; implicit none
#     include 'control_simulation.h'      
      integer :: mode
       
      select case(mode)

     case(1) ! initialization
!    Open positions file

        open(12,file="film_xmol",status="unknown",position="append")
        open(13,file="vel.dat",status="unknown",position="append")

    case(2)  ! writing of safe confs: conf_new and conf_xmol
!
! *** Restart configuration
!
! #### conf_new ####

         open(10,file="conf_new",status="unknown")
         write(10,'(i1)') 1
#ifdef STARS

         write(10,'(4i10)') n_mon,n_chain
         write(10, '(5i10)')n_mon_d,n_chain_d
         write(10,'(6i10)') n_mon_arm, n_arms, n_stars
#else

         write(10,'(4i10)') n_mon,n_chain,n_mon_d,n_chain_d
#endif

         write(10,120) n_cell_w_x,n_cell_w_y
         write(10,110) n_order
         do i_part = 1,n_part
             write(10,113)a_type(i_part),r0(1,i_part), r0(2,i_part),r0(3,i_part)
         end do
!obs         do i_order = 1,2
!obs             r_dummy = dt**i_order
             do i_part = 1,n_part
                 write(10,103) (v(i_dim,i_part),i_dim=1,n_dim)
             end do
             do i_part = 1,n_part
                 write(10,103) (a(i_dim,i_part),i_dim=1,n_dim)
             end do
!obs         end do
       !  if(f_explicit_wall) then !only if considering explicit wall atoms
#       if WALL == 1        
             do i_wall = 1,n_wall
                 write(10,103) (r_wall_equi(i_dim,i_wall),i_dim=1,n_dim)
             end do
       !  end if 
#       endif         
         write(10,103) (r0_twall(i_dim),i_dim=1,n_dim)
         do i_order = 1,n_order
             r_dummy = dt**i_order
             write(10,103) (rx_twall(i_dim,i_order)/r_dummy,i_dim=1,n_dim)
         end do
         write(10,130) (pbc_twall(i_dim),i_dim=1,n_dim)
         write(10,103) (r0_bwall(i_dim),i_dim=1,n_dim)
         write(10,130) (pbc_bwall(i_dim),i_dim=1,n_dim)
         write(10,103) (r0_spring_twall(i_dim),i_dim=1,n_dim)
         write(10,130) (s_force_grad(i_dim),i_dim=1,n_dim)
         close(10)

! #### conf_xmol ####

         open(11,file="conf_xmol",status="unknown")
         write(11,'(i5)') n_part
         write(11,*)
!    writes "brushes"      
         do i_part = 1,part_init_d
             if(a_type(i_part).eq.1) write(11,203) "Cl   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
             if(a_type(i_part).eq.2) write(11,203) "O   ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
         end do
!    writes droplet/melt
#ifndef PARTICLE_4
         do i_part = part_init_d+1,n_mon_tot
            write(11,203) "He   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
#else         
         do i_part = part_init_d+1,part_init_e
             write(11,203) "He   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
         do i_part = part_init_e+1,n_mon_tot
             write(11,203) "N    ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
#endif /*  particle 4 */        

! -----  Writes explicit wall atoms 

#       if WALL == 1 
             do i_part = n_mon_tot+1,n_part
                 write(11,203) " N  ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
             end do
#       endif
         close(11)
!         
! ----- Positions (xyz format) and velocities files 
!   
         case(3) 

! #### film_xmol ####

         write(12,'(i5)') n_part
         write(12,*)  " " !   i_time=  ",i_time
!brush
         do i_part = 1,part_init_d
#ifdef BRUSH_SEP
             if (i_part.le.n_chain*n_mon/2) then

             if(a_type(i_part).eq.1) write(12,203) "Cl   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)      !ori
             if(a_type(i_part).eq.2) write(12,203) "O   ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
             else

             if(a_type(i_part).eq.1) write(12,203) "Cl   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)      !ori
             if(a_type(i_part).eq.2) write(12,203) "Cu   ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
             endif
#else
             if(a_type(i_part).eq.1) write(12,203) "Cl   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)      !ori
             if(a_type(i_part).eq.2) write(12,203) "O   ",r0(1,i_part), r0(2,i_part),r0(3,i_part)

#endif
         end do
!melt         
#ifndef PARTICLE_4
#   if WALL != 1
         do i_part = part_init_d+1,n_mon_tot
#   else
         do i_part = part_init_d+1,n_mon_tot
#   endif
            write(12,203) "He   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
         
#else         
         do i_part = part_init_d+1,part_init_e
            write(12,203) "He   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
         do i_part = part_init_e+1,part_init_e+n_mon_e*n_chain_e
            write(12,203) "N    ",r0(1,i_part),r0(2,i_part),r0(3,i_part)
         end do
#ifdef STARS
          do i_part= part_init_e+n_mon_e*n_chain_e+1, n_mon_tot
             if(a_type(i_part).eq.5) write(12,203) "Al   ",r0(1,i_part),r0(2,i_part),r0(3,i_part)      !ori
             if(a_type(i_part).eq.6) write(12,203) "Pb   ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
         enddo
#endif
#endif /* defined particle 4 */        
!debug do i_part=1,n_mon_tot
!debug print*,a_type(i_part),r0(:,i_part)
!debug end do

!         if(f_explicit_wall) then
#       if WALL == 1 
             do i_part = n_mon_tot+1,n_part
                 write(12,203) " N  ",r0(1,i_part), r0(2,i_part),r0(3,i_part)
             end do
#       endif
!         end if

! #### vel.dat ####

         write(13,*) r_time
         write(13,*) "  "
         do i_part =1,n_mon_tot
              write(13,'(3f16.8)') v(:,i_part)
         end do

! -------------------------------
         case(4)    ! Store unfolded coordinates 

! #### film_xmol ####

         write(12,'(i5)') n_part
         write(12,*) " "  !,i_time
!brush
         do i_part = 1,part_init_d
             if(a_type(i_part).eq.1) write(12,203) "Cl   ",r0_unfold(1,i_part),r0_unfold(2,i_part),r0_unfold(3,i_part)      !ori
             if(a_type(i_part).eq.2) write(12,203) "O   ",r0_unfold(1,i_part), r0_unfold(2,i_part),r0_unfold(3,i_part)
         end do
!melt         
#ifndef PARTICLE_4
         do i_part = part_init_d+1,n_mon_tot
            write(12,203) "He   ",r0_unfold(1,i_part),r0_unfold(2,i_part),r0_unfold(3,i_part)
         end do
#else         
         do i_part = part_init_d+1,part_init_e
            write(12,203) "He   ",r0_unfold(1,i_part),r0_unfold(2,i_part),r0_unfold(3,i_part)
         end do
         do i_part = part_init_e+1,n_mon_tot
            write(12,203) "N    ",r0_unfold(1,i_part),r0_unfold(2,i_part),r0_unfold(3,i_part)
         end do
#endif /* defined particle 4 */        


!         if(f_explicit_wall) then
#       if WALL == 1
             do i_part = n_mon_tot+1,n_part
                 write(12,203) " N  ",r0_unfold(1,i_part), r0_unfold(2,i_part),r0_unfold(3,i_part)
             end do
#       endif

! #### vel.dat ####

         write(13,*) r_time
         write(13,*) "  "
         do i_part =1,n_mon_tot
              write(13,'(3f16.8)') v(:,i_part)
         end do
!      
      case(5) ! ------------ Closing of all the files and finish
          
#if STORE == 1 
!
! Drop the very last unfold_coordinate
!
      print '(/a/)', "* Writing last configuration to conf_unfold"
      open(unit=17,file='conf_unfold',status='unknown')
      do i_part = 1 , n_part 
      write(17,'(3e16.8)') r0_unfold(:,i_part)
      end do
      close(17)
#endif
      close(12) ! film_xmol
      close(13) !  vel.dat

 end select 
!claudio
#ifdef STARS

  103 format(3f25.7)
#else
  103 format(3f15.7)
#endif
  110 format(i3)
  113 format(i3,3f13.7)
  120 format(2i10)
  130 format(3i13)
  !203 format(a,3f13.4)
  203 format(a,3f13.9)
!
      end
