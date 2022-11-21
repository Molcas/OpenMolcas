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

subroutine PotIntd(Zeta,ZInv,rKappa,P,nZeta,la,lb,A,RB,Array,nArr,rFinal,lFinal,nDAO,CCoor,pot,ngrid,ncmp,DAO,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of potential (nordop=0)   *
!         or electric field (nordop > 0) on a grid                     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nArr, lFinal, nDAO, ngrid, ncmp, nOrdOp
real(kind=wp), intent(in) :: Zeta(nZeta), ZInv(nZeta), rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,*), DAO(*)
real(kind=wp), intent(inout) :: Array(nZeta*nArr), pot(ncmp,ngrid)
real(kind=wp), intent(out) :: rFinal(lFinal)
integer(kind=iwp) :: i, iAnga(4), icmp, igeo, ip1, ipc, ipIn, kab, lab, labcd, lcd, mabMax, mabMin, mArr, mcdMax, mcdMin, nComp, &
                     nFLOP, nMem, nT
real(kind=wp) :: CoorAC(3,2), Coori(3,4), TC(3)
logical(kind=iwp) :: NoSpecial
real(kind=wp), external :: ddot_
logical(kind=iwp), external :: EQ
external :: Fake, TNAI, XCff2D, XRys2D

nComp = nTri_Elem1(nOrdOp)

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
Coori(:,1) = A
Coori(:,2) = RB
mabMin = nTri3_Elem1(max(la,lb)-1)
mabMax = nTri3_Elem1(la+lb)-1
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
mcdMin = nTri3_Elem1(nOrdOp-1)
mcdMax = nTri3_Elem1(nOrdOp)-1
lab = (mabMax-mabMin+1)
kab = nTri_Elem1(la)*nTri_Elem1(lb)
lcd = (mcdMax-mcdMin+1)
labcd = lab*lcd

! Compute FLOP's and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Distribute the work array

ip1 = 1+nZeta*max(labcd,lcd*nMem)
mArr = nArr-max(labcd,lcd*nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if

nT = nZeta
NoSpecial = .true.

do igeo=1,ngrid
  do i=1,3
    TC(i) = CCoor(i,igeo)
    CoorAC(i,2) = TC(i)
    Coori(i,3) = TC(i)
    Coori(i,4) = TC(i)
  end do

  ! Compute integrals with the Rys quadrature.

  call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coori,CoorAC,mabMin,mabMax,mcdMin,mcdMax, &
           Array(ip1),mArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

  if (nOrdOp == 0) then

    ! Use the HRR to compute the required primitive integrals.

    call HRR(la,lb,A,RB,Array(ip1),nZeta,nMem,ipIn)
    ipc = ip1-1+ipin

    ! Trace primitive integrals with density

    pot(1,igeo) = pot(1,igeo)-ddot_(nDAO,Array(ipc),1,DAO,1)
  else
    call DGetMO(Array(ip1),nZeta*lab,nZeta*lab,lcd,Array,lcd)
    call HRR(la,lb,A,RB,Array,lcd*nZeta,nMem,ipIn)
    call DGetMO(Array(ipIn),lcd,lcd,nZeta*kab,rFinal,nZeta*kab)
    ipc = 1
    do icmp=1,ncomp
      pot(icmp,igeo) = pot(icmp,igeo)+ddot_(nDAO,rFinal(ipc),1,DAO,1)
      ipc = ipc+nDAO
    end do
  end if
end do

return

end subroutine PotIntd
