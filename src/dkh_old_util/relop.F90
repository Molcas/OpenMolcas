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

implicit real*8(A-H,O-Z)
#include "crelop.fh"

!write(6,100)
!100 format(/,' ****** RELATIVISTIC OPERATORS V 1.0 - BERND HESS ******'//)
PI = 4.d0*atan(1.d0)
ZWP = 2.d0*PI
ZWPH32 = ZWP**1.5d0
ZWPH12 = sqrt(ZWP)
SQPI = sqrt(PI)
VELIT = CLightAU
PREA = 1.d0/(VELIT*VELIT)
CSQ = VELIT*VELIT
FAK(1) = 1.d0
!GHALB(1) = SQPI
do I=2,26
  !GHALB(I) = GHALB(I-1)*(dble(I)-1.5D0)
  FAK(I) = FAK(I-1)*dble(I-1)
end do

! BINOMIALKOEFFIZIENTEN

IMAX = 20
BCO(1) = 1.d0
IBIAS = 1
JBIAS = 1
K = IMAX-1
do I=1,K
  ADD = 0.d0
  do J=1,I
    JBIAS = JBIAS+1
    BCO(JBIAS) = ADD+BCO(IBIAS)
    ADD = BCO(IBIAS)
    IBIAS = IBIAS+1
  end do
  JBIAS = JBIAS+1
  BCO(JBIAS) = 1.d0
end do

do N=1,20
  GA(N) = GAM(N-1)
end do

return

end subroutine RELOP
