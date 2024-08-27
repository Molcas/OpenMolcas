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
function EstI(Zeta,rKapAB,nAlpha,nBeta,Coeff1,niBas,Coeff2,njBas,xab,nab,Scrt,nScrt,IndZ)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Constants, only: Zero

implicit none
real*8 EstI
integer nAlpha, nBeta, niBas, njBas, nab, nScrt
real*8 Zeta(nAlpha*nBeta), rKapAB(nAlpha,nBeta), Coeff1(nAlpha,niBas), Coeff2(nBeta,njBas), xab(nAlpha*nBeta), Scrt(nScrt)
integer IndZ(nAlpha*nBeta+1)
integer mZeta, iZeta, iAlpha, iBeta, iEta, iDelta, iGamma, iBas, jBas, ijBas, iHigh
integer, external :: iDAMax_
real*8 rab, rcd, AInt

#ifdef _DEBUGPRINT_
write(6,*) 'Esti:mZeta=',IndZ(nAlpha*nBeta)
call RecPrt('Esti:xab',' ',xab,1,nAlpha*nBeta)
call RecPrt('Esti:Coeff1',' ',Coeff1,nAlpha,niBas)
call RecPrt('Esti:Coeff2',' ',Coeff2,nBeta,njBas)
#endif

mZeta = IndZ(nAlpha*nBeta+1)
call dcopy_(niBas*njBas,[Zero],0,Scrt,1)
do iZeta=1,mZeta
  iBeta = (IndZ(iZeta)-1)/nAlpha+1
  iAlpha = IndZ(iZeta)-(iBeta-1)*nAlpha
  rab = xab(iZeta)
  do iEta=1,mZeta
    iDelta = (IndZ(iEta)-1)/nAlpha+1
    iGamma = IndZ(iEta)-(iDelta-1)*nAlpha
    rcd = xab(iEta)
    AInt = rab*rcd
    do iBas=1,niBas
      do jBas=1,njBas
        ijBas = (jBas-1)*niBas+iBas
        Scrt(ijBas) = Scrt(ijBas)+abs(Coeff1(iAlpha,iBas)*Coeff2(iBeta,jBas))*abs(Coeff1(iGamma,iBas)*Coeff2(iDelta,jBas))*AInt
      end do
    end do
  end do
end do
iHigh = iDAMax_(niBas*njBas,Scrt,1)
EstI = sqrt(Scrt(iHigh))

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Zeta)
  call Unused_real_array(rKapAB)
  call Unused_integer(nab)
end if

end function EstI
