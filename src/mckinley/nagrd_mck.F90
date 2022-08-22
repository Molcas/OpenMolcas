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

use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp), parameter :: nHess = 1, nPAO = 1 ! Hess, PAO: Dummy arrays
integer(kind=iwp) :: iAlpha, iAnga(4), iBeta, iCnt, iDCRT(0:7), iElem, iIrrep, indi, Indx(3,4), ipA, ipAOff, ipB, ipBOff, &
                     iuvwx(4), iZeta, JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), kCnt, kCnttp, kdc, kndgrd(3,4,0:7), lDCRT, LmbdT, &
                     mOp(4), nArray, nb, nDCRT, nGr, nip, nRys
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), Fact, Hess(nHess), PAO(nPAO), TC(3), tfac
logical(kind=iwp) :: DiffCnt, jfg(4), JfGrd(3,4), JfHss(4,3,4,3), kfgrd(3,4), Tr(4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ
external :: Cff2D, Fake, TNAI1
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!iRout = 150
!iPrint = nPrint(iRout)

!if (iPrint >= 99) then
!  write(u6,*) ' In NAGrd: nArr=',nArr
!end if

nRys = nHer

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
if (nip-1 > nArr) write(u6,*) ' nip-1 > nArr'
nArray = nArr-nip+1

iIrrep = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
! Dummies
call ICopy(144*nIrrep,[0],0,JndHss,1)
call LCopy(144,[.false.],0,jfHss,1)

call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)
if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
iuvwx(1) = iu
iuvwx(2) = iv
mOp(1) = nOp(1)
mOp(2) = nOp(2)

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

! Loop over nuclear centers

nb = nZeta*nElem(la)*nElem(lb)
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
    call LCopy(12,[.false.],0,JFgrd,1)
    call ICopy(12*nIrrep,[0],0,jndGrd,1)
    do iCnt=1,2
      JfGrd(iDCar,iCnt) = IfGrad(iDCar,iCnt)
    end do
    do ICnt=1,2
      if (IfGrad(idcar,iCnt)) then
        do iIrrep=0,nIrrep-1
          jndGrd(iDCar,iCnt,iIrrep) = IndGrd(iIrrep)
        end do
      end if
    end do

    Tr(1) = .false.
    Tr(2) = .false.
    Tr(3) = .false.
    Tr(4) = .false.
    if ((kdc+kCnt) == iDCnt) then
      Tr(3) = .true.
      JfGrd(iDCar,1) = .true.
      JfGrd(iDCar,2) = .true.
      do iIrrep=0,nIrrep-1
        jndGrd(iDCar,3,iIrrep) = -IndGrd(iIrrep)
      end do
    end if

    do lDCRT=0,nDCRT-1
      call lCopy(12,JfGrd,1,kfGrd,1)
      call iCopy(12*nIrrep,JndGrd,1,kndgrd,1)
      mOp(3) = NrOpr(iDCRT(lDCRT))
      mOp(4) = mOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coora(1,3),1)
      call dcopy_(3,TC,1,Coora(1,4),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      if (Eq(A,RB) .and. EQ(A,TC)) cycle
      if (EQ(A,TC)) then
        kfGrd(iDCar,1) = .false.
        do iIrrep=0,nIrrep-1
          kndGrd(iDCar,1,iirrep) = 0
        end do
      end if
      if (EQ(RB,TC)) then
        kfGrd(iDCar,2) = .false.
        do iIrrep=0,nIrrep-1
          kndgrd(iDCar,2,iIrrep) = 0
        end do
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

      do iElem=1,nElem(la)*nElem(lb)*ngr
        do iZeta=1,nZeta
          tfac = Two*rKappa(iZeta)*Pi*ZInv(iZeta)
          indi = (iElem-1)*nZeta+iZeta
          Array(nip+indi-1) = tfac*Array(nip+indi-1)
        end do
      end do

#     ifdef _DEBUGPRINT_
      call RecPrt('In NaGrd PI',' ',Array(nip),nb,3)
      call RecPrt('In NaGrd PI',' ',rFinal,nb,nrOp)
#     endif
      call SmAdNa(Array(nip),nb,rFinal,mop,loper,KndGrd,iuvwx,kfGrd,Indx,idcar,Fact,JFG,tr)
      !if (iPrint > 23) call RecPrt('In NaGrd FI',' ',rFinal,nb,nrOp)
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_logical_array(Trans)
end if

end subroutine NAGrd_mck
