    module mpi_module
contains

subroutine mpi_initialize()
#ifdef MPI    
    use mpi_commons
    call MPI_INIT(ierr)
#endif
end subroutine mpi_initialize
!!! not yet !  ----------------------------------------------------------    
!!! not yet         subroutine initial_broadcast()
!!! not yet #       include 'control_simulation.h'
!!! not yet 
!!! not yet ! Broadcast to all other procs what 0 has learned (I/O)
!!! not yet 
!!! not yet         use commons
!!! not yet         use ziggurat ! gaussian number generator suite
!!! not yet         use mpi_commons
!!! not yet         implicit none
!!! not yet         logical, parameter :: debug=.true.
!!! not yet 
!!! not yet #ifdef MPI
!!! not yet 
!!! not yet ! Times
!!! not yet        call MPI_BCAST(n_relax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_obser,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet ! Dimensions 
!!! not yet        call MPI_BCAST(n_mon_tot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_part,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_loop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) !WARN nloop could be removed
!!! not yet        call MPI_BCAST(n_chain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_chain_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_mon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_mon_d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet ! Some flags
!!! not yet        !call MPI_BCAST(solv_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(interact_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet !       call MPI_BCAST(n_dim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet ! force switch-on
!!! not yet        call MPI_BCAST(r_2_min_time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_2_min_init,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(boundary,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet ! Binning variables
!!! not yet !print*,my_rank,boundary,r_2_min_time,r_2_min_init
!!! not yet 
!!! not yet        call MPI_BCAST(n_neigh_fl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_bin_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_bin_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_bin_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_cell,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(n_max_cell,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_bin_x,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_bin_y,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_bin_z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet ! Thermostats, cutoff, others
!!! not yet        call MPI_BCAST(friction,n_type,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(skin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(skin_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(z_space_wall,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(sig,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(w_d,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(w_r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet ! DEBUG
!!! not yet !print*,my_rank,skin,skin_2
!!! not yet 
!!! not yet ! Interaction variables
!!! not yet !       LJ       
!!! not yet        call MPI_BCAST(range_2,n_type**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(e_shift,n_type**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(sigma_2,n_type**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(epsil,n_type**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(mass_type,n_type,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet !  Bonded interactions
!!! not yet        call MPI_BCAST(k_chain,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_chain_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(r_chain,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet ! Wall interactions 
!!! not yet        call MPI_BCAST(a_w,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(sigma_w,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet ! Temperature control
!!! not yet        call MPI_BCAST(sig,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(temp_final,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet      ! DPD 
!!! not yet        call MPI_BCAST(w_d,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(w_r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(mz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet        call MPI_BCAST(kptr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet 
!!! not yet ! External force
!!! not yet        call MPI_BCAST(const_force,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet 
!!! not yet ! Global variables
!!! not yet 
!!! not yet 
!!! not yet        allocate (g_force(3,n_mon_tot)) 
!!! not yet        g_force(:,:) = 0 
!!! not yet       
!!! not yet 
!!! not yet ! Allocate the variables in all processors but 0  
!!! not yet 
!!! not yet 
!!! not yet         if (my_rank /= 0 ) then
!!! not yet 
!!! not yet !**** Dynamic allocation of variables  for non 0 procs. ****
!!! not yet 
!!! not yet ! Position and forces:
!!! not yet 
!!! not yet         allocate (   r0(n_dim,n_part),   r0_old(n_dim,n_part), spabs(n_dim,n_part), force(n_dim,n_part)  ) 
!!! not yet         allocate( r_head_old(n_chain*n_mon,3) ) 
!!! not yet         allocate ( f_on_heads(3,n_chain) )
!!! not yet ! veloc., accel. and higher derivatives
!!! not yet        allocate(rx(n_part,n_dim,n_order))  
!!! not yet ! Mass and atom types, vectors:
!!! not yet 	allocate(   a_type(n_part),              &
!!! not yet                     mass(n_part),                &
!!! not yet                     random(n_part)   )         ! random vector
!!! not yet ! binning
!!! not yet 	allocate  ( mic_count(n_dim,n_part),mic_old(n_dim,n_part) )
!!! not yet 	allocate (  ff_list(0:n_neigh_fl,n_mon_tot) )
!!! not yet ! DPD forces
!!! not yet         allocate (force_d(3,n_mon_tot) , force_r(3,n_mon_tot) )
!!! not yet 
!!! not yet ! Now all the vectors should be allocated.
!!! not yet 
!!! not yet !---- Initialize random number generators in the rest of the procs
!!! not yet !
!!! not yet       call inr250(mz,iseed,kptr)
!!! not yet 
!!! not yet !      
!!! not yet ! --- Initialize random number generator used for random force in DPD
!!! not yet !          (ziggurat)
!!! not yet 
!!! not yet       call zigset(iseed)
!!! not yet 
!!! not yet 	
!!! not yet    end if  ! proc 0
!!! not yet 
!!! not yet ! Broadcasting of the current values of the vectors 
!!! not yet 
!!! not yet          call MPI_bcast(r0(:,:)    ,n_part*n_dim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet          call MPI_bcast(r0_old(:,:),n_part*n_dim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet          call MPI_bcast(spabs(:,:) ,n_part*n_dim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet          call MPI_bcast(f_on_heads ,n_chain*n_dim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet !         call MPI_bcast(force(:,:) ,n_part*n_dim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet          force(:,:) =0
!!! not yet 
!!! not yet          call MPI_bcast(r_head_old,n_chain*n_mon*3, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!! not yet 
!!! not yet          call MPI_bcast(dt_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!! not yet          call MPI_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!! not yet 
!!! not yet          call MPI_bcast(rx,n_part*n_dim*n_order,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
!!! not yet 
!!! not yet          call MPI_bcast (a_type,n_part,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet          call MPI_bcast (mass,n_part,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
!!! not yet          call MPI_bcast (random,n_part,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet          call MPI_bcast(mic_count,n_part*n_dim,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet          call MPI_bcast(mic_old,n_part*n_dim,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!! not yet 
!!! not yet          call MPI_bcast(ff_list,n_mon_tot*(1+n_neigh_fl),MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
!!! not yet 
!!! not yet #endif
!!! not yet ! MPI
!!! not yet end subroutine initial_broadcast
!!! not yet 
!!! not yet 
!!! not yet    
!!! not yet            subroutine interval(kk,n1,n2)
!!! not yet !kk = my_rank      
!!! not yet         use commons
!!! not yet         use mpi_commons
!!! not yet       
!!! not yet         implicit none
!!! not yet       
!!! not yet       
!!! not yet         integer,intent(in) ::  kk
!!! not yet         integer ,intent(out):: n1,n2
!!! not yet         integer npart,nsup
!!! not yet       
!!! not yet         if(mod(n_mon_tot,nproc).eq.0) then
!!! not yet            n1=1+kk*n_mon_tot/nproc
!!! not yet            n2=(kk+1)*n_mon_tot/nproc
!!! not yet         else
!!! not yet            if(kk.lt.nproc-1) then
!!! not yet               n1=1+kk*int(n_mon_tot/nproc)
!!! not yet               n2=(kk+1)*int(n_mon_tot/nproc)
!!! not yet            else
!!! not yet               n1=1+kk*int(n_mon_tot/nproc)
!!! not yet               n2=(kk+1)*int(n_mon_tot/nproc)+mod(n_mon_tot,nproc)
!!! not yet            endif
!!! not yet         endif
!!! not yet       
!!! not yet    return
!!! not yet    
!!! not yet    end subroutine interval
!!! not yet 
!!! not yet !ori CEM    
!!! not yet !ori CEM    subroutine share(x)
!!! not yet !ori CEM    
!!! not yet !ori CEM      use commons, only : NATOM,ndim,i,j,k
!!! not yet !ori CEM      use mpicommons
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      implicit none
!!! not yet !ori CEM    
!!! not yet !ori CEM      include 'parameters.inc'
!!! not yet !ori CEM    
!!! not yet !ori CEM      real x(NATOM,ndim)
!!! not yet !ori CEM    
!!! not yet !ori CEM      real transferx(NATOM*ndim)
!!! not yet !ori CEM    
!!! not yet !ori CEM      integer n1,n2
!!! not yet !ori CEM      integer nc(0:nproc-1)
!!! not yet !ori CEM      integer ndis(0:nproc-1)
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      do i=0,nproc-1
!!! not yet !ori CEM         call interval(i,n1,n2)
!!! not yet !ori CEM         nc(i)=(n2-n1+1)*3
!!! not yet !ori CEM         ndis(i)=(n1-1)*3
!!! not yet !ori CEM      enddo
!!! not yet !ori CEM    
!!! not yet !ori CEM      call interval(my_rank,n1,n2)
!!! not yet !ori CEM    
!!! not yet !ori CEM      k=0
!!! not yet !ori CEM      do i=n1,n2
!!! not yet !ori CEM         do j=1,ndim
!!! not yet !ori CEM            k=k+1
!!! not yet !ori CEM            transferx(k)=x(i,j)
!!! not yet !ori CEM         enddo
!!! not yet !ori CEM      enddo
!!! not yet !ori CEM     
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      call
!!! not yet !ori CEM    MPI_GATHERV(transferx,nc,MPI_DOUBLE_PRECISION,transferx,nc,ndis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      if(my_rank.eq.0) then
!!! not yet !ori CEM         k=0
!!! not yet !ori CEM         do i=1,NATOM
!!! not yet !ori CEM            do j=1,ndim
!!! not yet !ori CEM               k=k+1
!!! not yet !ori CEM               x(i,j)=transferx(k)
!!! not yet !ori CEM            enddo
!!! not yet !ori CEM         enddo
!!! not yet !ori CEM      endif
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      call MPI_BCAST(x,ndim*NATOM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!!! not yet !ori CEM    
!!! not yet !ori CEM    
!!! not yet !ori CEM      return
!!! not yet !ori CEM    
!!! not yet !ori CEM    end subroutine share
!!! not yet 
!!! not yet 
end module mpi_module
