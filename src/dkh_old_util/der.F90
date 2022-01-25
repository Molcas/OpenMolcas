!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

real*8 function DER(IDER,IS1,IS2,AL,BE)
! CALCULATE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS

implicit real*8(A-H,O-Z)
integer IS1(3), IS2(3), I1(2,3), I2(2,3)
real*8 F(2), G(2)
#include "crelop.fh"

do I=1,3
  do J=1,2
    I1(J,I) = IS1(I)
    I2(J,I) = IS2(I)
  end do
end do
I1(1,IDER) = I1(1,IDER)+1
I1(2,IDER) = I1(2,IDER)-1
I2(1,IDER) = I2(1,IDER)+1
I2(2,IDER) = I2(2,IDER)-1
L1 = IS1(IDER)+1
goto(10,11,12,13,14),L1
101 write(6,100) IDER,IS1,IS2,AL,BE
100 format(' ILLEGAL ANGULAR MOMENTUM (DER)'/,' IDER,IS1,IS2,AL,BE PRINTED'/,1X,7I5,3X,2d20.8)
call Abend()
10 F(1) = -2.d0*AL
J1 = 1
goto 19

11 F(2) = 1.d0
F(1) = -2.d0*AL
J1 = 2
goto 19

12 F(2) = 2.d0
F(1) = -2.d0*AL
J1 = 2
goto 19

13 F(2) = 3.d0
F(1) = -2.d0*AL
J1 = 2
goto 19

14 F(2) = 4.0d0
F(1) = -2.0d0*AL
J1 = 2

19 L2 = IS2(IDER)+1
goto(20,21,22,23,24),L2
goto 101

20 G(1) = -2.d0*BE
J2 = 1
goto 29

21 G(2) = 1.d0
G(1) = -2.d0*BE
J2 = 2
goto 29

22 G(2) = 2.d0
G(1) = -2.d0*BE
J2 = 2
goto 29

23 G(2) = 3.d0
G(1) = -2.d0*BE
J2 = 2
goto 29

24 G(2) = 4.0d0
G(1) = -2.0d0*BE
J2 = 2

29 SUM = 0.d0
do I=1,J1
  do J=1,J2
    II = I1(I,1)+I2(J,1)
    JJ = I1(I,2)+I2(J,2)
    KK = I1(I,3)+I2(J,3)
    ANG = THETA(II+JJ,KK)*PHI(JJ,II)
    if (ANG == 0.d0) goto 31
    EX = -dble(II+JJ+KK+2)*0.5d0
    SUM = SUM+F(I)*G(J)*0.5d0*ANG*GA(II+JJ+KK+2)*(AL+BE)**EX
31  continue
  end do
end do
DER = SUM

return

end function DER
