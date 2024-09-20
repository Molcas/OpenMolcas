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
subroutine EFPrm(Zeta,ZInv,rKappa,P,rFinal,nZeta,nComp,la,lb,A,RB,Array,nArr,Ccoor,nOrdOp)
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

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, rys2d_kernel, tval_kernel
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nZeta, nComp, la, lb, nArr, nOrdOp
real(kind=wp), intent(in) :: Zeta(nZeta), ZInv(nZeta), rKappa(nZeta), P(nZeta,3), A(3), RB(3), Ccoor(3)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp)
real(kind=wp), intent(inout) :: Array(nZeta*nArr)
integer(kind=iwp) :: iAnga(4), ip1, ip2, ip3, ipIn, kab, lab, labcd, lcd, mabMax, mabMin, mArr, mcdMax, mcdMin, nFlop, nMem, nT
real(kind=wp) :: CoorAC(3,2), Coori(3,4)
logical(kind=iwp) :: NoSpecial
procedure(cff2d_kernel) :: XCff2D
procedure(modu2_kernel) :: Fake
procedure(rys2d_kernel) :: XRys2D
procedure(tval_kernel) :: TNAI
logical(kind=iwp), external :: EQ
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iElem, jElem
character(len=80) :: Label
#endif

rFinal(:,:,:,:) = Zero

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
Coori(:,1) = A(:)
Coori(:,2) = RB(:)
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
mcdMin = nTri3_Elem1(nOrdOp-1)
mcdMax = nTri3_Elem1(nOrdOp)-1
lab = (mabMax-mabMin+1)
kab = nTri_Elem1(la)*nTri_Elem1(lb)
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
  CoorAC(:,1) = A(:)
else
  CoorAC(:,1) = RB(:)
end if

CoorAC(:,2) = Ccoor(:)
Coori(:,3) = Ccoor(:)
Coori(:,4) = Ccoor(:)

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

call DGetMO(Array(ip3),lcd,lcd,nZeta*kab,rFinal,nZeta*kab)
rFinal(:,:,:,1:lcd) = -rFinal(:,:,:,1:lcd)

#ifdef _DEBUGPRINT_
write(u6,*) ' In EFPrm la,lb=',la,lb
do iElem=1,nTri_Elem1(la)
  do jElem=1,nTri_Elem1(lb)
    if (lcd == 1) then
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: rFinal (',iElem,',',jElem,') '
      call RecPrt(Label,' ',rFinal(1,iElem,jElem,1),nZeta,1)
    else if (lcd == 3) then
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: rFinal (',iElem,',',jElem,',x) '
      call RecPrt(Label,' ',rFinal(1,iElem,jElem,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: rFinal (',iElem,',',jElem,',y) '
      call RecPrt(Label,' ',rFinal(1,iElem,jElem,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' EFPrm: rFinal (',iElem,',',jElem,',z) '
      call RecPrt(Label,' ',rFinal(1,iElem,jElem,3),nZeta,1)
    end if
  end do
end do
#endif

return

end subroutine EFPrm
