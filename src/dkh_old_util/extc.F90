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

function EXTC(LAMBDA,AL,BE,L1,M1,N1,L2,M2,N2)

use crelop, only: GA, IMAX
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: EXTC
integer(kind=iwp), intent(in) :: LAMBDA, L1, M1, N1, L2, M2, N2
real(kind=wp), intent(in) :: AL, BE
integer(kind=iwp) :: II, IS1(3), IS2(3), JJ, KK
real(kind=wp) :: ANG, EX, OV1, OV2, SUM_
real(kind=wp), external :: DER, PHI, THETA

! CALCULATE ANGULAR AND RADIAL PART

II = L1+L2
JJ = M1+M2
KK = N1+N2
IMAX = II+JJ+KK+3
if (IMAX > 20) then

  ! ERROR BRANCH: ANGULAR MOMENTUM  > MAXIMUM GIVEN BY ARRAY GA

  write(u6,1002) L1,M1,N1,L2,M2,N2,LAMBDA
  call Abend()

end if

! COMPUTE INTEGRAL OVER DERIVATIVE OF THE FUNCTIONS

IS1(1) = L1
IS1(2) = M1
IS1(3) = N1
IS2(1) = L2
IS2(2) = M2
IS2(3) = N2
SUM_ = DER(1,IS1,IS2,AL,BE)+DER(2,IS1,IS2,AL,BE)+DER(3,IS1,IS2,AL,BE)

! NORMALIZATION

II = L1+L1
JJ = M1+M1
KK = N1+N1
ANG = THETA(II+JJ,KK)*PHI(JJ,II)
EX = -Half*real(II+JJ+KK+3,kind=wp)
OV1 = Half*ANG*GA(II+JJ+KK+3)*((AL+AL)**EX)
!write(u6,*) ' EXTC OV1',L1,M1,N1,AL,ANG,sqrt(1/OV1)
II = L2+L2
JJ = M2+M2
KK = N2+N2
ANG = THETA(II+JJ,KK)*PHI(JJ,II)
EX = -Half*real(II+JJ+KK+3,kind=wp)
OV2 = Half*ANG*GA(II+JJ+KK+3)*((BE+BE)**EX)
!write(u6,*) ' EXTC OV1',L2,M2,N2,BE,ANG,sqrt(1/OV2)
EXTC = SUM_/sqrt(OV1*OV2)

return

1002 format(' ILLEGAL ANGULAR MOMENTUM (PVP)'/,' L1,M1,N1,L2,M2,N2,LAMBDA PRINTED'/,1X,7I5)

end function EXTC
