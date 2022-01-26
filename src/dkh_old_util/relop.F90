!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Bernd Artur Hess                                 *
!***********************************************************************

subroutine RELOP()
! SUBROUTINE RELOP INITIALIZES THE COMMON BLOCK USED BY
! THE RELOP PACKAGE
! V 1.0 - 12.3.86 - BERND HESS

use DKH_Info, only: CLightAU
use Constants, only: Zero, One, Two, OneHalf, Pi
use Definitions, only: wp, iwp

implicit none
#include "crelop.fh"
integer(kind=iwp) :: I, IBIAS, J, JBIAS, K, N
real(kind=wp) :: ADD
real(kind=wp), external :: GAM

!write(u6,100)
!100 format(/,' ****** RELATIVISTIC OPERATORS V 1.0 - BERND HESS ******'//)
ZWP = Two*PI
ZWPH32 = ZWP**OneHalf
ZWPH12 = sqrt(ZWP)
SQPI = sqrt(PI)
VELIT = CLightAU
PREA = One/(VELIT*VELIT)
CSQ = VELIT*VELIT
FAK(1) = One
!GHALB(1) = SQPI
do I=2,26
  !GHALB(I) = GHALB(I-1)*(real(I,kind=wp)-OneHalf)
  FAK(I) = FAK(I-1)*real(I-1,kind=wp)
end do

! BINOMIALKOEFFIZIENTEN

IMAX = 20
BCO(1) = One
IBIAS = 1
JBIAS = 1
K = IMAX-1
do I=1,K
  ADD = Zero
  do J=1,I
    JBIAS = JBIAS+1
    BCO(JBIAS) = ADD+BCO(IBIAS)
    ADD = BCO(IBIAS)
    IBIAS = IBIAS+1
  end do
  JBIAS = JBIAS+1
  BCO(JBIAS) = One
end do

do N=1,20
  GA(N) = GAM(N-1)
end do

return

end subroutine RELOP
