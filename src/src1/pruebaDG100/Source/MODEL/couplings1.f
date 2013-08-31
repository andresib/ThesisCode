ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER  (PI=3.141592653589793D0)


      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'

      GC_1 = -(EE*COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*EE*COMPLEXI)/3.000000D+00
      GC_16 = (CKM22*EE*COMPLEXI)/(SW*SQRT__2)
      GC_21 = -(CW*EE*COMPLEXI)/(2.000000D+00*SW)
      GC_22 = (CW*EE*COMPLEXI)/(2.000000D+00*SW)
      GC_23 = -(EE*COMPLEXI*SW)/(6.000000D+00*CW)
      END
