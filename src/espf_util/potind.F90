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

subroutine PotIntd(Zeta,ZInv,rKappa,P,nZeta,la,lb,A,RB,Array,nArr,final,lFinal,nDAO,CCoor,pot,ngrid,ncmp,DAO,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of potential (nordop=0)   *
!         or electric field (nordop > 0) on a grid                     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!***********************************************************************

use k2_arrays
implicit real*8(A-H,O-Z)
! Used for normal nuclear attraction integrals
external TNAI, Fake, XCff2D, XRys2D
#include "real.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
!#include "print.fh"
real*8 Zeta(nZeta), ZInv(nZeta), rKappa(nZeta), P(nZeta,3), A(3), RB(3), CCoor(3,*), Array(nZeta*nArr), pot(ncmp,ngrid), DAO(*), &
       final(lfinal)
! Local arrays
!real*8 C(3), TC(3), Coori(3,4), CoorAC(3,2)
real*8 TC(3), Coori(3,4), CoorAC(3,2)
logical EQ, NoSpecial
integer iAnga(4)
!character ChOper(0:7)*3
!data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement functions
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

nComp = (nOrdOp+1)*(nOrdOp+2)/2

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

! Compute FLOP's and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Distribute the work array

ip1 = 1+nZeta*max(labcd,lcd*nMem)
mArr = nArr-max(labcd,lcd*nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
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
    call DGetMO(Array(ipIn),lcd,lcd,nZeta*kab,final,nZeta*kab)
    ipc = 1
    do icmp=1,ncomp
      pot(icmp,igeo) = pot(icmp,igeo)+ddot_(nDAO,final(ipc),1,DAO,1)
      ipc = ipc+nDAO
    end do
  end if
end do

return

end subroutine PotIntd
