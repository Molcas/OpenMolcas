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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine RdMx(RadMax,Exp,nExp,Cff,nCff,cdMax,EtMax)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             August '91                                               *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 exp(nExp), Cff(nExp,nCff)
#include "real.fh"
#ifdef _DEBUGPRINT_
#include "print.fh"
#endif

#ifdef _DEBUGPRINT_
iRout = 201
iPrint = nPrint(iRout)

call RecPrt('Exp',' ',Exp,nExp,1)
call RecPrt('Cff',' ',Cff,nExp,nCff)
#endif
do iExp=1,nExp

  cc = DDot_(nCff,Cff(iExp,1),nExp,Cff(iExp,1),nExp)
  c = sqrt(cc)

  Alpha = exp(iExp)
  Beta = exp(iExp)
  Zeta = Alpha+Beta
  if (Zeta > Zero) then
    Eta = Alpha+Beta
    Rho = (Zeta*Eta)/(Zeta+Eta)

    ssss = c**4*Two*sqrt(Rho/Pi)*(Pi/Zeta)**(Three/Two)*(Pi/Eta)**(Three/Two)
    if (sqrt(ssss) > RadMax) then
      RadMax = sqrt(ssss)
      EtMax = Eta
      cdMax = sqrt(ssss)
    end if
  end if

end do

end subroutine RdMx
