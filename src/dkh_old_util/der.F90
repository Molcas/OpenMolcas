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

function DER(IDER,IS1,IS2,AL,BE)
! CALCULATE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS

use crelop, only: GA
use Constants, only: Zero, One, Two, Three, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: DER
integer(kind=iwp), intent(in) :: IDER, IS1(3), IS2(3)
real(kind=wp), intent(in) :: AL, BE
integer(kind=iwp) :: I, I1(2,3), I2(2,3), II, J, J1, J2, JJ, KK, L1, L2
real(kind=wp) :: ANG, EX, F(2), G(2), SUM_
real(kind=wp), external :: PHI, THETA

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
J1 = 0
select case (L1)
  case default
    write(u6,100) IDER,IS1,IS2,AL,BE
    call Abend()
  case (1)
    F(1) = -TWo*AL
    J1 = 1
  case (2)
    F(2) = One
    F(1) = -Two*AL
    J1 = 2
  case (3)
    F(2) = TWo
    F(1) = -TWo*AL
    J1 = 2
  case (4)
    F(2) = Three
    F(1) = -Two*AL
    J1 = 2
  case (5)
    F(2) = Four
    F(1) = -Two*AL
    J1 = 2
end select

L2 = IS2(IDER)+1
J2 = 0
select case (L2)
  case default
    write(u6,100) IDER,IS1,IS2,AL,BE
    call Abend()
  case (1)
    G(1) = -Two*BE
    J2 = 1
  case (2)
    G(2) = One
    G(1) = -Two*BE
    J2 = 2
  case (3)
    G(2) = Two
    G(1) = -Two*BE
    J2 = 2
  case (4)
    G(2) = Three
    G(1) = -Two*BE
    J2 = 2
  case (5)
    G(2) = Four
    G(1) = -Two*BE
    J2 = 2
end select

SUM_ = Zero
do I=1,J1
  do J=1,J2
    II = I1(I,1)+I2(J,1)
    JJ = I1(I,2)+I2(J,2)
    KK = I1(I,3)+I2(J,3)
    ANG = THETA(II+JJ,KK)*PHI(JJ,II)
    if (ANG == Zero) cycle
    EX = -real(II+JJ+KK+2,kind=wp)*Half
    SUM_ = SUM_+F(I)*G(J)*Half*ANG*GA(II+JJ+KK+2)*(AL+BE)**EX
  end do
end do
DER = SUM_

return

100 format(' ILLEGAL ANGULAR MOMENTUM (DER)'/,' IDER,IS1,IS2,AL,BE PRINTED'/,1X,7I5,3X,2d20.8)

end function DER
