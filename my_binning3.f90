
subroutine my_binning3()
  
#include 'control_simulation.h'
      use commons ; implicit none
      logical,parameter :: binning_debug=.false.
  
  
    ff_list(:,1:n_part) = 0 
    binpnt(:) = 0
    bin(:) = 0

    do i_part=1,n_part !n_mon_tot

!     if(pntcopy(i).eq.0) cycle

        ix = int ( r0(1,i_part) / r_bin_x  )  
        iy = int ( r0(2,i_part) / r_bin_y  ) 
        iz = int ( r0(3,i_part) / r_bin_z  ) 

        ib=iz*NCELLY*NCELLX+iy*NCELLX+ix+1

        bin(i)=binpnt(ib)
        binpnt(ib) = i
     

!!!     if(ib.lt.1.or.ib.gt.NCELL) then
!!!        write(6,*) my_rank
!!!        write(6,*) i,ix,iy,iz
!!!        write(6,*) rlocp(i,1),rlocp(i,2),rlocp(i,3)
!!!!        stop
!!!     endif

  enddo
  
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        do i_part=1,n_part
!  do i=1,nlocal
    

!ii=ii+1

!     if(ii.gt.NTOT) cycle

!     xtmp=rloc(i,1)
!     ytmp=rloc(i,2)
!     ztmp=rloc(i,3)

        r0_tmp(:) = r0(:,i_part)
     
!     ixx=(xtmp+0.5*DISTX)/CX
!     iyy=(ytmp+0.5*DISTY)/CY
!     if(qwall) then
!        izz=ztmp/CZ
!     else
!        izz=(ztmp+0.5*DISTZ)/CZ
!     endif
        ixx = int ( r0(1,i_part) / r_bin_x )  
        iyy = int ( r0(2,i_part) / r_bin_y ) 
        izz = int ( r0(3,i_part) / r_bin_z ) 

     do k = 0,26
         ix = ixx + mod(k,3) - 1
         if (ix.lt.0) ix = NCELLX- 1
         if (ix.eq.NCELLX) ix = 0
         iy = iyy + mod(k/3,3) - 1
         if (iy.lt.0) iy = NCELLY- 1
         if (iy.eq.NCELLY) iy = 0
         iz = izz + k/9 - 1
         if (iz.lt.0) iz = NCELLZ- 1
         if (iz.eq.NCELLZ) iz = 0
         ib = iz*NCELLY*NCELLX+iy*NCELLX+ ix + 1
         j_part = binpnt(ib)

         j_part = binpnt(ib)

10          if(j_part /= 0) then

!             jj=pntcopy(j)

             ! CPU P_ii takes all the atoms
             if(mod(my_rank,psquare+1).ne.0) then


                 if(whereami(my_rank)) then

                     ! CPU P_ij takes the odd combinations

                     if(mod(ii+jj,2).eq.0) then
                         goto 35
                     endif

                 else
                     ! CPU P_ji takes the even combinations

                     if(mod(ii+jj,2).ne.0) then
                         goto 35
                     endif

                 endif

             else

                 if(j.le.i) then
                     goto 35
                 endif


             endif
           
! don't compute substrat-substrat neighbours
           
           if(ntype(ii).eq.2.and.jj.ne.0) then
              if(ntype(jj).eq.2) then
                 goto 35
              endif
           endif

           dr(1)=xtmp-rlocp(j,1)
           dr(2)=ytmp-rlocp(j,2)
           dr(3)=ztmp-rlocp(j,3)
           
           
           if(dr(1).gt.0.5*DISTX) then
              dr(1)=dr(1)-DISTX
           elseif(dr(1).lt.-0.5*DISTX) then
              dr(1)=dr(1)+DISTX
           endif

           if(dr(2).gt.0.5*DISTY) then
              dr(2)=dr(2)-DISTY
           elseif(dr(2).lt.-0.5*DISTY) then
              dr(2)=dr(2)+DISTY
           endif
           
           if(qperz) then
              if(dr(3).gt.0.5*DISTZ) then
                 dr(3)=dr(3)-DISTZ
              elseif(dr(3).lt.-0.5*DISTZ) then
                 dr(3)=dr(3)+DISTZ
              endif
           endif
           
           d=0.
           do l=1,ndim
              d=d+dr(l)*dr(l)
           enddo
           d=SQRT(d)
           
           if(d.lt.rv.and.ii.ne.jj) then
              list(i,0)=list(i,0)+1
              list(i,list(i,0))=j
           endif
           
35         j = bin(j)

           goto 10

        endif        
     enddo
  enddo
  


  t2=MPI_WTime()
  tneigh=tneigh+t2-t1
  

  return
  
end subroutine my_binning3
