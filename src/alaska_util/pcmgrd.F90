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
! Copyright (C) 1995,2001, Roland Lindh                                *
!***********************************************************************

subroutine PCMgrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, May '95                                 *
!                                                                      *
!             Modified to PCM gradients September 2001, Lund, by       *
!             R. Lindh.                                                *
!***********************************************************************

use PCM_arrays, only: PCM_SQ, PCMTess
use Center_Info, only: dc
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iAnga(4), iBeta, iCar, iDAO, iDCRT(0:7), ipA, ipAOff, ipB, ipBOff, ipDAO, iPrint, iRout, &
                     iStb(0:7), iTs, iuvwx(4), iZeta, j, JndGrd(3,4), lDCRT, LmbdT, lOp(4), mGrad, mRys, nArray, nDAO, nDCRT, &
                     nDiff, nip, nStb
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), EInv, Eta, Fact, Q, TC(3)
logical(kind=iwp) :: JfGrad(3,4)
!character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: NrOpr
external :: Fake, TNAI1, XCff2D
#include "print.fh"
#include "rctfld.fh"
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

#include "macros.fh"
unused_var(rFinal)
unused_var(nHer)
unused_var(Ccoor(1))
unused_var(nComp)

iRout = 151
iPrint = nPrint(iRout)

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nElem(la)*nElem(lb)*nElem(nOrdOp)
if (nip-1 > nZeta*nArr) then
  write(u6,*) 'nip-1 > nZeta*nArr'
  call Abend()
end if
nArray = nZeta*nArr-nip+1

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = nOrdOp
iAnga(4) = 0
Coori(:,1) = A(:)
Coori(:,2) = RB(:)
if (la >= lb) then
  CoorAC(:,1) = A(:)
else
  CoorAC(:,1) = RB(:)
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)

ipAOff = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
  ipBOff = ipBOff+1
end do

! Modify the density matrix with the prefactor

nDAO = nElem(la)*nElem(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

! Generate stabilizer of C, i.e. a center of a tile.

nStb = 1
iStb(0) = 0

! Loop over the tiles

do iTs=1,nTs
  Q = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
  if (Q == Zero) cycle
  ! Pick up the tile coordinates
  C(1:3) = PCMTess(1:3,iTs)

  if (iPrint >= 99) call RecPrt('C',' ',C,1,3)
  call DCR(LmbdT,iStabM,nStabM,iStb,nStb,iDCRT,nDCRT)
  Fact = -Q*real(nStabM,kind=wp)/real(LmbdT,kind=wp)

  Array(ipDAO:ipDAO+nZeta*nDAO-1) = Fact*pack(DAO,.true.)

  iuvwx(3) = nStb
  iuvwx(4) = nStb
  JndGrd(:,1:2) = IndGrd(:,:)
  do i=1,3
    do j=1,2
      JfGrad(i,j) = IfGrad(i,j)
    end do
  end do

  ! No derivatives with respect to the third or fourth center.
  ! The positions of the points in the external field are frozen.

  JndGrd(:,3) = 0
  JfGrad(1,3) = .false.
  JfGrad(2,3) = .false.
  JfGrad(3,3) = .false.
  JndGrd(:,4) = 0
  JfGrad(1,4) = .false.
  JfGrad(2,4) = .false.
  JfGrad(3,4) = .false.
  mGrad = 0
  do iCar=1,3
    do i=1,2
      if (JfGrad(iCar,i)) mGrad = mGrad+1
    end do
  end do
  if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
  if (mGrad == 0) cycle

  do lDCRT=0,nDCRT-1
    lOp(3) = NrOpr(iDCRT(lDCRT))
    lOp(4) = lOp(3)
    call OA(iDCRT(lDCRT),C,TC)
    CoorAC(:,2) = TC(:)
    Coori(:,3) = TC(:)
    Coori(:,4) = TC(:)

    ! Compute integrals with the Rys quadrature.

    nDiff = 1
    mRys = (la+lb+2+nDiff+nOrdOp)/2
    Eta = One
    EInv = One
    call Rysg1(iAnga,mRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[Eta],[EInv],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
               Array(nip),nArray,TNAI1,Fake,XCff2D,Array(ipDAO),nDAO*nElem(nOrdOp),Grad,nGrad,JfGrad,JndGrd,lOp,iuvwx)

    !call RecPrt(' In PCMgrd:Grad',' ',Grad,nGrad,1)
  end do  ! End loop over DCRs

end do    ! End loop over centers in the external field

return

end subroutine PCMgrd
