! This version is made for using with the r250
! packages. Initialization of them must also be made.
! from the calling program.


!     ALGORITHM 488 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN COMM. ACM, VOL. 17, NO. 12,
!     P. 704.                                                           
      FUNCTION GRAND(N)                                                 
! EXCEPT ON THE FIRST CALL GRAND RETURNS A
! PSEUDO-RANDOM NUMBER HAVING A GAUSSIAN (I.E.
! NORMAL) DISTRIBUTION WITH ZERO MEAN AND UNIT
! STANDARD DEVIATION.  THUS, THE DENSITY IS  F(X) =
! EXP(-0.5*X**2)/SQRT(2.0*PI). THE FIRST CALL
! INITIALIZES GRAND AND RETURNS ZERO.
! THE PARAMETER N IS DUMMY.
! GRAND CALLS A FUNCTION RAND, AND IT IS ASSUMED THAT
! SUCCESSIVE CALLS TO RAND(0) GIVE INDEPENDENT
! PSEUDO- RANDOM NUMBERS DISTRIBUTED UNIFORMLY ON (0,
! 1), POSSIBLY INCLUDING 0 (BUT NOT 1).
! THE METHOD USED WAS SUGGESTED BY VON NEUMANN, AND
! IMPROVED BY FORSYTHE, AHRENS, DIETER AND BRENT.
! ON THE AVERAGE THERE ARE 1.37746 CALLS OF RAND FOR
! EACH CALL OF GRAND.
! WARNING - DIMENSION AND DATA STATEMENTS BELOW ARE
!           MACHINE-DEPENDENT.
! DIMENSION OF D MUST BE AT LEAST THE NUMBER OF BITS
! IN THE FRACTION OF A FLOATING-POINT NUMBER.
! THUS, ON MOST MACHINES THE DATA STATEMENT BELOW
! CAN BE TRUNCATED.
! IF THE INTEGRAL OF SQRT(2.0/PI)*EXP(-0.5*X**2) FROM
! A(I) TO INFINITY IS 2**(-I), THEN D(I) = A(I) -
! A(I-1).  
      use commons, only: mz,kptr
      real :: rand(2)
      DIMENSION D(60)
      DATA D(1), D(2), D(3), D(4), D(5), D(6), D(7),                    &
     & D(8), D(9), D(10), D(11), D(12), D(13),                          &
     & D(14), D(15), D(16), D(17), D(18), D(19),                        &
     & D(20), D(21), D(22), D(23), D(24), D(25),                        &
     & D(26), D(27), D(28), D(29), D(30), D(31),                        &
     & D(32) /0.674489750,0.475859630,0.383771164,                      &
     & 0.328611323,0.291142827,0.263684322,                             &
     & 0.242508452,0.225567444,0.211634166,                             &
     & 0.199924267,0.189910758,0.181225181,                             &
     & 0.173601400,0.166841909,0.160796729,                             &
     & 0.155349717,0.150409384,0.145902577,                             &
     & 0.141770033,0.137963174,0.134441762,                             &
     & 0.131172150,0.128125965,0.125279090,                             &
     & 0.122610883,0.120103560,0.117741707,                             &
     & 0.115511892,0.113402349,0.111402720,                             &
     & 0.109503852,0.107697617/
      DATA D(33), D(34), D(35), D(36), D(37), D(38),                    &
     & D(39), D(40), D(41), D(42), D(43), D(44),                        &
     & D(45), D(46), D(47), D(48), D(49), D(50),                        &
     & D(51), D(52), D(53), D(54), D(55), D(56),                        &
     & D(57), D(58), D(59), D(60)                                       &
     & /0.105976772,0.104334841,0.102766012,                            &
     & 0.101265052,0.099827234,0.098448282,                             &
     & 0.097124309,0.095851778,0.094627461,                             &
     & 0.093448407,0.092311909,0.091215482,                             &
     & 0.090156838,0.089133867,0.088144619,                             &
     & 0.087187293,0.086260215,0.085361834,                             &
     & 0.084490706,0.083645487,0.082824924,                             &
     & 0.082027847,0.081253162,0.080499844,                             &
     & 0.079766932,0.079053527,0.078358781,                             &
     & 0.077681899/
! END OF MACHINE-DEPENDENT STATEMENTS
! U MUST BE PRESERVED BETWEEN CALLS.
      DATA U /0.0/
! INITIALIZE DISPLACEMENT A AND COUNTER I.
      A = 0.0
      I = 0
! INCREMENT COUNTER AND DISPLACEMENT IF LEADING BIT
! OF U IS ONE.
   10 U = U + U
      IF (U.LT.1.0) GO TO 20
      U = U - 1.0
      I = I + 1
      A = A - D(I)
      GO TO 10
! FORM W UNIFORM ON 0 .LE. W .LT. D(I+1) FROM U.
   20 W = D(I+1)*U
! FORM V = 0.5*((W-A)**2 - A**2). NOTE THAT 0 .LE. V
! .LT. LOG(2).
      V = W*(0.5*W-A)
! GENERATE NEW UNIFORM U.
        call r250(mz,rand,2,2,kptr)
   30 U = rand(1) ! ori RAND(0)
! ACCEPT W AS A RANDOM SAMPLE IF V .LE. U.
      IF (V.LE.U) GO TO 40
! GENERATE RANDOM V.
      V = rand(2) ! ori RAND(0)
! LOOP IF U .GT. V.
      IF (U.GT.V) GO TO 30
! REJECT W AND FORM A NEW UNIFORM U FROM V AND U.
      U = (V-U)/(1.0-U)
      GO TO 20
! FORM NEW U (TO BE USED ON NEXT CALL) FROM U AND V.
   40 U = (U-V)/(1.0-V)
! USE FIRST BIT OF U FOR SIGN, RETURN NORMAL VARIATE.
      U = U + U
      IF (U.LT.1.0) GO TO 50
      U = U - 1.0
      GRAND = W - A
      RETURN
   50 GRAND = A - W
      RETURN
        END FUNCTION GRAND 

