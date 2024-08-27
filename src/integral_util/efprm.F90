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
! Copyright (C) 1991,1995, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine EFPrm(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,final,nZeta,nComp,la,lb,A,RB,nRys,Array,nArr,Ccoor,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of electric field         *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
! Modified for explicit code, R. Lindh, February '95.                  *
!***********************************************************************

use Constants, only: Zero, One

implicit none
external TNAI, Fake, XCff2D, XRys2D
integer nZeta, la, lb, nComp, nAlpha, nBeta, nRys, nArr, nOrdOp
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta), rKappa(nZeta), &
       P(nZeta,3), A(3), RB(3), Array(nZeta*nArr), Ccoor(3)
real*8 Coori(3,4), CoorAC(3,2)
logical EQ, NoSpecial
integer iAnga(4)
integer ixyz, nElem, nabSz
integer mabMin, mabMax, mcdMin, mcdMax, lab, kab, lcd, labcd, ip1, ip2, nMem, mArr, nT, ip3, ipIn, nFlop
#ifdef _DEBUGPRINT_
integer iElem, jElem
character*80 Label
#endif
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

#ifdef _DEBUGPRINT_
call RecPrt(' In EFPrm: Alpha',' ',Alpha,nAlpha,1)
call RecPrt(' In EFPrm: Beta',' ',Beta,nBeta,1)
#else
if (.false.) call Unused_Real_Array(Alpha)
if (.false.) call Unused_Real_Array(Beta)
#endif

final(:,:,:,:) = Zero

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)
mabMin = nabSz(max(la,lb)-1)+1
mabMax = nabSz(la+lb)
if (EQ(A,RB)) mabMin = nabSz(la+lb-1)+1
mcdMin = nabSz(nOrdOp-1)+1
mcdMax = nabSz(nOrdOp)
lab = (mabMax-mabMin+1)
kab = nElem(la)*nElem(lb)
lcd = (mcdMax-mcdMin+1)
labcd = lab*lcd

! Compute Flop's and size of work array which HRR will Use.

call mHRR(la,lb,nFLOP,nMem)

! Distribute the work array

ip2 = 1
ip1 = ip2+nZeta*max(labcd,lcd*nMem)
mArr = nArr-max(labcd,lcd*nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if

call dcopy_(3,Ccoor,1,CoorAC(1,2),1)
call dcopy_(3,Ccoor,1,Coori(1,3),1)
call dcopy_(3,Ccoor,1,Coori(1,4),1)

! Compute integrals with the Rys-Gauss quadrature.

nT = nZeta
NoSpecial = .true.
call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,CCoor,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
         Array(ip1),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

! The integrals are now ordered as ijkl,e,f

! a) Change the order to f,ijkl,e
! b) Unfold e to ab, f,ijkl,ab
! c) Change the order back to ijkl,ab,f

!a)--

call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array(ip2),lcd)

!b)-- Use the HRR to unfold e to ab

call HRR(la,lb,A,RB,Array(ip2),lcd*nZeta,nMem,ipIn)
ip3 = ip2-1+ipIn

!c)--

call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,final,nZeta*kab)
call DScal_(nZeta*kab*lcd,-One,final,1)

#ifdef _DEBUGPRINT_
write(6,*) ' In EFPrm la,lb=',la,lb
do iElem=1,nElem(la)
  do jElem=1,nElem(lb)
    if (lcd == 1) then
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: Final (',iElem,',',jElem,') '
      call RecPrt(Label,' ',final(1,iElem,jElem,1),nZeta,1)
    else if (lcd == 3) then
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: Final (',iElem,',',jElem,',x) '
      call RecPrt(Label,' ',final(1,iElem,jElem,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: Final (',iElem,',',jElem,',y) '
      call RecPrt(Label,' ',final(1,iElem,jElem,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: Final (',iElem,',',jElem,',z) '
      call RecPrt(Label,' ',final(1,iElem,jElem,3),nZeta,1)
    end if
  end do
end do
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nRys)

end subroutine EFPrm
