! Random numbers generator 
      SUBROUTINE R250(MZ,RAN,n_dim,N,KPTR)
!
! STORES N REAL RANDOM NUMBERS IN RAN USING THE INTEGERS IN MZ
! THE RANDOM NUMBERS HAVE UNIFORM DISTRIBUTION IN THE INTERVAL (0,1)
!
      PARAMETER(RINV=1./2147483647.)
      INTEGER MZ(250)
      INTEGER N,KPTR,n_dim
      REAL RAN(n_dim)
      INTEGER I,L
!
!
      I = 0
      L = N + KPTR - 1
!
 10   CONTINUE
!
      KMIN = KPTR
      KMAX = MIN(L,147)
      DO 100 K = KMIN,KMAX
        MZ(K)  = IEOR(MZ(K),MZ(K + 103))
        I      = I + 1
        RAN(I) = MZ(K) * RINV
 100  CONTINUE
!
      KMIN = MAX(KMIN,KMAX+1)
      KMAX = MIN(L,250)
!VOCL LOOP,NOVREC
      DO 200 K = KMIN,KMAX
        MZ(K)  = IEOR(MZ(K),MZ(K - 147))
        I      = I + 1
        RAN(I) = MZ(K) * RINV
 200  CONTINUE
!
      IF  (KMAX .EQ. 250) THEN
        KPTR = 1
        L    = L - 250
        GOTO 10
      END IF
!
      KPTR = KMAX + 1
!
      RETURN
      END
