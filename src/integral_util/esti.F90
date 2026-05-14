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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
function EstI(nAlpha,nBeta,Coeff1,niBas,Coeff2,njBas,xab,Scrt,nScrt,IndZ)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp) :: EstI
integer(kind=iwp), intent(in) :: nAlpha, nBeta, niBas, njBas, nScrt, IndZ(nAlpha*nBeta+1)
real(kind=wp), intent(in) :: Coeff1(nAlpha,niBas), Coeff2(nBeta,njBas), xab(nAlpha*nBeta)
real(kind=wp), intent(out) :: Scrt(nScrt)
integer(kind=iwp) :: iAlpha, iBas, iBeta, iDelta, iEta, iGamma, iHigh, ijBas, iZeta, jBas, mZeta
real(kind=wp) :: A_Int, rab, rcd
integer(kind=iwp), external :: iDAMax_

#ifdef _DEBUGPRINT_
write(u6,*) 'Esti:mZeta=',IndZ(nAlpha*nBeta)
call RecPrt('Esti:xab',' ',xab,1,nAlpha*nBeta)
call RecPrt('Esti:Coeff1',' ',Coeff1,nAlpha,niBas)
call RecPrt('Esti:Coeff2',' ',Coeff2,nBeta,njBas)
#endif

mZeta = IndZ(nAlpha*nBeta+1)
Scrt(1:niBas*njBas) = Zero
do iZeta=1,mZeta
  iBeta = (IndZ(iZeta)-1)/nAlpha+1
  iAlpha = IndZ(iZeta)-(iBeta-1)*nAlpha
  rab = xab(iZeta)
  do iEta=1,mZeta
    iDelta = (IndZ(iEta)-1)/nAlpha+1
    iGamma = IndZ(iEta)-(iDelta-1)*nAlpha
    rcd = xab(iEta)
    A_Int = rab*rcd
    do iBas=1,niBas
      do jBas=1,njBas
        ijBas = (jBas-1)*niBas+iBas
        Scrt(ijBas) = Scrt(ijBas)+abs(Coeff1(iAlpha,iBas)*Coeff2(iBeta,jBas))*abs(Coeff1(iGamma,iBas)*Coeff2(iDelta,jBas))*A_Int
      end do
    end do
  end do
end do
iHigh = iDAMax_(niBas*njBas,Scrt,1)
EstI = sqrt(Scrt(iHigh))

return

end function EstI
