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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine NAGrd_mck( &
#                    define _CALLING_
#                    include "grd_mck_interface.fh"
                    )
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the nuclear attraction integrals. *
!          Something is wrong here                                     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             October 1991                                             *
!              Anders Bernhardsson 1995                                *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp), parameter :: nHess = 1, nPAO = 1 ! Hess, PAO: Dummy arrays
integer(kind=iwp) :: iAnga(4), iBeta, iCnt, iDCRT(0:7), iElem, indi, Indx(3,4), ipA, ipAOff, ipB, ipBOff, iuvwx(4), &
                     JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), kCnt, kCnttp, kdc, kndgrd(3,4,0:7), lDCRT, LmbdT, mOp(4), nArray, nb, &
                     nDCRT, nGr, nip, nRys
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), Fact, Hess(nHess), PAO(nPAO), TC(3)
logical(kind=iwp) :: DiffCnt, jfg(4), JfGrd(3,4), JfHss(4,3,4,3), kfgrd(3,4), Tr(4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
external :: Cff2D, Fake, TNAI1

#include "macros.fh"
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(Trans)

!iRout = 150
!iPrint = nPrint(iRout)

!if (iPrint >= 99) then
!  write(u6,*) ' In NAGrd_McK: nArr=',nArr
!end if

nRys = nHer

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
if (nip-1 > nArr) write(u6,*) ' nip-1 > nArr'
nArray = nArr-nip+1

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
! Dummies
JndHss(:,:,:,:,0:nIrrep-1) = 0
JfHss(:,:,:,:) = .false.

Coora(:,1) = A
Coora(:,2) = RB
Coori(:,1) = A
Coori(:,2) = RB
if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if
iuvwx(1) = iu
iuvwx(2) = iv
mOp(1) = nOp(1)
mOp(2) = nOp(2)

ipAOff = ipA
do iBeta=1,nBeta
  Array(ipAOff:ipAOff+nAlpha-1) = Alpha
  ipAOff = ipAOff+nAlpha
end do

ipBOff = ipB
do iBeta=1,nBeta
  Array(ipBOff:ipBOff+nAlpha-1) = Beta(iBeta)
  ipBOff = ipBOff+nAlpha
end do

! Loop over nuclear centers

nb = nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (kCnttp == iCnttp_Dummy) cycle
  if (dbsc(kCnttp)%Charge == Zero) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
    DiffCnt = (IfGrad(iDCar,1) .or. IfGrad(iDCar,2))
    if ((.not. DiffCnt) .and. (kdc+kCnt /= iDCnt)) cycle

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    !Fact = -dbsc(kCnttp)%Charge*real(nStabM*nIrrep,kind=wp)/real(LmbdT*dc(kdc+kCnt)%nStab,kind=wp)
    Fact = -dbsc(kCnttp)%Charge*real(nStabM,kind=wp)/real(LmbdT,kind=wp)
    !if (iPrint >= 99) then
    !  write(u6,*) ' Charge=',dbsc(kCnttp)%Charge
    !  Write(u6,*) 'NZeta=',nzeta
    !  write(u6,*) 'NrOp=',nrop
    !  write(u6,*) ' Fact=',Fact
    !end if
    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab
    JndGrd(:,:,0:nIrrep-1) = 0
    JfGrd(:,:) = .false.
    JfGrd(iDCar,1:2) = IfGrad(iDCar,1:2)
    do ICnt=1,2
      if (IfGrad(idcar,iCnt)) JndGrd(iDCar,iCnt,0:nIrrep-1) = IndGrd(0:nIrrep-1)
    end do

    Tr(1) = .false.
    Tr(2) = .false.
    Tr(3) = .false.
    Tr(4) = .false.
    if ((kdc+kCnt) == iDCnt) then
      Tr(3) = .true.
      JfGrd(iDCar,1:2) = .true.
      JndGrd(iDCar,3,0:nIrrep-1) = -IndGrd(0:nIrrep-1)
    end if

    do lDCRT=0,nDCRT-1
      kndgrd(:,:,0:nIrrep-1) = JndGrd(:,:,0:nIrrep-1)
      kfgrd(:,:) = JfGrd(:,:)
      mOp(3) = NrOpr(iDCRT(lDCRT))
      mOp(4) = mOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      CoorAC(:,2) = TC
      Coora(:,3) = TC
      Coora(:,4) = TC
      Coori(:,3) = TC
      Coori(:,4) = TC
      if (Eq(A,RB) .and. EQ(A,TC)) cycle
      if (EQ(A,TC)) then
        kfGrd(iDCar,1) = .false.
        kndGrd(iDCar,1,0:nIrrep-1) = 0
      end if
      if (EQ(RB,TC)) then
        kfGrd(iDCar,2) = .false.
        kndgrd(iDCar,2,0:nIrrep-1) = 0
      end if

      if (kfGrd(idcar,1)) then
        JFG(1) = .true.
      else
        JFG(1) = .false.
      end if
      if (kfGrd(idcar,2)) then
        JFG(2) = .true.
      else
        JFG(2) = .false.
      end if
      JFG(3) = .false.
      JFG(4) = .false.
      call Rysg2(iAnga,nRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coora,CoorAC, &
                 Array(nip),nArray,TNAI1,Fake,Cff2D,PAO,nPAO,Hess,nHess,kfGrd,kndGrd,JfHss,JndHss,mOp,iuvwx,Jfg,nGr,Indx,.true., &
                 .false.,tr)

      do iElem=1,nTri_Elem1(la)*nTri_Elem1(lb)*ngr
        indi = nip+(iElem-1)*nZeta
        Array(indi:indi+nZeta-1) = Two*rKappa(:)*Pi*ZInv(:)*Array(indi:indi+nZeta-1)
      end do

#     ifdef _DEBUGPRINT_
      call RecPrt('In NaGrd_McK PI',' ',Array(nip),nb,3)
      call RecPrt('In NaGrd_McK PI',' ',rFinal,nb,nrOp)
#     endif
      call SmAdNa(Array(nip),nb,rFinal,mop,loper,KndGrd,iuvwx,Indx,idcar,Fact,tr)
      !if (iPrint > 23) call RecPrt('In NaGrd_McK FI',' ',rFinal,nb,nrOp)
    end do
  end do
end do

return

end subroutine NAGrd_mck
