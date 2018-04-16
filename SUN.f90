!
! ---------------------------------------------------------------------
!
      SUBROUTINE SUN(ISEED,RANF,N)
!
! STORES N REAL RANDOM NUMBERS IN RANF
! ISEED IS A STARTUP VALUE
!
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!
      INTEGER ISEED
      REAL RANF(N)
!
!
!
      DATA FACTOR /41475557.0D0/, TWO28 /268435456.0D0/
!
!
!
      IF (ISEED .GE. 0) THEN
        R=DBLE(ISEED)/TWO28
        R=DMOD(R*FACTOR,1.0D0)
        ISEED=-1
      END IF
!
       DO 100 I = 1,N
        R=DMOD(R*FACTOR,1.0D0)
        RANF(I) = SNGL(R)
 100  CONTINUE
!
      RETURN
      END
