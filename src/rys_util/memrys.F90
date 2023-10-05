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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine MemRys(iAnga,MemPrm)
! This routine will compute the memory requirement of RYS
! Memory requirement is per primitive!

use Gateway_global, only: FMM_shortrange
use Index_Functions, only: nTri3_Elem1
use Breit, only: nOrdOp, nComp
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iAnga(4)
integer(kind=iwp), intent(out) :: MemPrm
integer(kind=iwp) :: la, labcd, labMax, labMin, lb, lB00, lB01, lB10, lc, lcdMax, lcdMin, ld, nabcd, nabMax, ncdMax, nRys, nabcdN

la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)

nRys = (la+lb+lc+ld+2)/2
if (nOrdOp == 0) then
  nRys = (la+lb+lc+ld+2)/2
  !nRys = (la+lb+lc+ld+4)/2
else if (nOrdOp == 1) then
  nRys = (la+lb+lc+ld+4)/2
else if (nOrdOp == 2) then
  nRys = (la+lb+lc+ld+4)/2
end if

labMin = nTri3_Elem1(max(la,lb)-1)
labMax = nTri3_Elem1(la+lb)-1
lcdMin = nTri3_Elem1(max(lc,ld)-1)
lcdMax = nTri3_Elem1(lc+ld)-1
labcd = (labMax-labMin+1)*(lcdMax-lcdMin+1)
#ifdef _DEBUGPRINT_
write(u6,*) ' labMin=',labMin
write(u6,*) ' labMax=',labMax
write(u6,*) ' lcdMin=',lcdMin
write(u6,*) ' lcdMax=',lcdMax
#endif
MemPrm = 0
! [a0|c0]
!  6 elements in the case of integrals for spin-spin coupling
MemPrm = MemPrm+nComp*labcd
!                                                                      *
!***********************************************************************
!                                                                      *
! For FMM, we only want short-range integrals, using twice the memory
! to store full and long-range components (which are subtracted)
! This option is not active for nOrdOp/=0

if (FMM_shortrange) MemPrm = MemPrm+nComp*labcd
!                                                                      *
!***********************************************************************
!                                                                      *
if (nOrdOp == 0) then
  nabMax = la+lb+nOrdOp
  !nabMin = max(la,lb)
  ncdMax = lc+ld+nOrdOp
  !ncdMin = max(lc,ld)
else
  nabMax = la+lb+2
  !nabMin = max(la,lb)
  ncdMax = lc+ld+2
  !ncdMin = max(lc,ld)
end if
nabcd = (nabMax+1)*(ncdMax+1)      ! ordinary 2D integrals
nabcdN = (nabMax-1+1)*(ncdMax-1+1)  ! extended 2D integrals
lB10 = max(min(nabMax-1,1),0)
lB01 = max(min(ncdMax-1,1),0)
lB00 = max(min(min(nabMax,ncdMax),1),0)
! Normalization
MemPrm = MemPrm+1
! Ordinary 2D-Integrals
MemPrm = MemPrm+nabcd*3*nRys
! Extended 2D-integrals
if (nOrdOp /= 0) MemPrm = MemPrm+nabcdN*3*2*nRys
! Coefficients for recurrence relations
MemPrm = MemPrm+3*nRys+3*nRys+3*nRys*(lB10+lB01+lB00)
! Roots
MemPrm = MemPrm+nRys
! The inverse of the arguments
MemPrm = MemPrm+1
! Arguments
MemPrm = MemPrm+1
! Expanded versions of Zeta, ZetInv, Eta, EtaInv, rKapab, rKapcd, P and Q
MemPrm = MemPrm+12
#ifdef _DEBUGPRINT_
write(u6,*) ' [e0|f0] integrals   :',labcd
write(u6,*) ' Normalization factor:',1
write(u6,*) ' 2D-integrals        :',nabcd*3*nRys
if (nOrdOp /= 0) write(u6,*) ' 2D-integrals extend :',nabcdN*3*nRys
write(u6,*) ' PAQP vector         :',3*nRys
write(u6,*) ' QCPQ vector         :',3*nRys
write(u6,*) ' B10 coefficients    :',nRys*3*lB10
write(u6,*) ' B00 coefficients    :',nRys*3*lB00
write(u6,*) ' B01 coefficients    :',nRys*3*lB01
write(u6,*) ' Roots               :',nRys
write(u6,*) ' Inverse arguments   :',1
write(u6,*) ' Arguments           :',1
#endif

return

end subroutine MemRys
