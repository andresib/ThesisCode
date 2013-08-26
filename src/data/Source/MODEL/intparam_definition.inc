ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        CKM22 = 1.000000D+00

        COMPLEXI = (0.000000D+00,1.000000D+00)

        MZ__EXP__2 = MZ**2

        MZ__EXP__4 = MZ**4

        SQRT__2 = DSQRT(2.000000D+00)

        MH__EXP__2 = MH**2

        AEW = 1.000000D+00/AEWM1

        MW = DSQRT(MZ__EXP__2/2.000000D+00+DSQRT(MZ__EXP__4/4.000000D
     $   +00-(AEW*PI*MZ__EXP__2)/(GF*SQRT__2)))

        SQRT__AEW = DSQRT(AEW)

        EE = 2.000000D+00*SQRT__AEW*DSQRT(PI)

        MW__EXP__2 = MW**2

        SW2 = 1.000000D+00-MW__EXP__2/MZ__EXP__2

        CW = DSQRT(1.000000D+00-SW2)

        SQRT__SW2 = DSQRT(SW2)

        SW = SQRT__SW2

        G1 = EE/CW

        GW = EE/SW

        V = (2.000000D+00*MW*SW)/EE

        V__EXP__2 = V**2

        LAM = MH__EXP__2/(2.000000D+00*V__EXP__2)

        YB = (YMB*SQRT__2)/V

        YT = (YMT*SQRT__2)/V

        YTAU = (YMTAU*SQRT__2)/V

        MUH = DSQRT(LAM*V__EXP__2)

        GW__EXP__2 = GW**2

        CW__EXP__2 = CW**2

        EE__EXP__2 = EE**2

        SW__EXP__2 = SW**2

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      AS = G**2/4/PI

      SQRT__AS = DSQRT(AS)

      G__EXP__2 = G**2

C     
C     Definition of the EW coupling used in the write out of aqed
C     
      GAL(1) = 1D0
      GAL(2) = 1D0

