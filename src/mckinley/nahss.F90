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
! Copyright (C) 1991, Anders Bernhardsson                              *
!               1991, Roland Lindh                                     *
!***********************************************************************

subroutine NAHss( &
#                define _CALLING_
#                include "hss_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the gradient of the nuclear attraction integrals. *
!                                                                      *
!     Author: Anders Bernhardsson & Roland Lindh,                      *
!             Dept. of Theoretical Chemistry, University               *
!             of Lund, SWEDEN.                                         *
!             October 1991                                             *
!***********************************************************************

use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp
use Center_Info, only: dc
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
#include "hss_interface.fh"
#include "Molcas.fh"
#include "disp.fh"
#include "disp2.fh"
integer(kind=iwp) :: iAlpha, iAnga(4), iAtom, iBeta, iCar, iCent, iComp, iDAO, iDCRT(0:7), iIrrep, Indx(3,4), ipA, ipAOff, ipArr, &
                     ipB, ipBOff, ipDAO, iStop, iuvwx(4), iZeta, jAtom, jCar, JndGrd(0:2,0:3,0:7), JndHss(0:3,0:2,0:3,0:2,0:7), &
                     kCar, kCent, kCnt, kCnttp, kdc, lDCRT, LmbdT, Maxi, Mini, mOp(4), nArray, nDAO, nDCRT, nDisp, nFinal, nip, &
                     nnIrrep, nRys
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, TC(3)
logical(kind=iwp) :: IfG(0:3), JfGrd(0:2,0:3), JfHss(0:3,0:2,0:3,0:2), Tr(0:3)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ, TF
external :: TNAI1, Fake, Cff2D
! Statement functions
integer(kind=iwp) :: nElem, ixyz, itri, i1, i2
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
itri(i1,i2) = max(i1,i2)*(max(i1,i2)-1)/2+min(i1,i2)

#ifdef _DEBUGPRINT_
write(u6,*) ' In NAHss: nArr=',nArr
#endif

nRys = nHer

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipDAO = nip
nip = nip+nAlpha*nBeta*nElem(la)*nElem(lb)
if (nip-1 > nArr) then
  write(u6,*) 'NAHss: nip-1 > nArr'
  write(u6,*) 'nip,nArr=',nip,nArr
  call Abend()
end if
ipArr = nip
nArray = nArr-nip+1

iIrrep = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
call dcopy_(3,A,1,Coori(1,1),1)
call dcopy_(3,RB,1,Coori(1,2),1)
if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
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

! Modify the density matrix with the prefactor

nDAO = nElem(la)*nElem(lb)
do iDAO=1,nDAO
  do iZeta=1,nZeta
    Fact = Two*rkappa(iZeta)*Pi/Zeta(iZeta)
    DAO(iZeta,iDAO) = Fact*DAO(iZeta,iDAO)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('DAO',' ',DAO,nZeta,nDAO)
#endif

! Loop over nuclear centers

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (kCnttp == iCnttp_Dummy) cycle
  if (dbsc(kCnttp)%Charge == Zero) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = -dbsc(kCnttp)%Charge*real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    call DYaX(nZeta*nDAO,Fact,DAO,1,Array(ipDAO),1)

    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    do lDCRT=0,nDCRT-1

      mOp(3) = NrOpr(iDCRT(lDCRT))
      mOp(4) = mOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      if (EQ(A,TC) .and. EQ(A,RB)) cycle

      ! Initialize JfGrd, JndGrd, JfHss, and JndHss.

      call LCopy(12,[.false.],0,JfGrd,1)
      call ICopy(nSym*4*3,[0],0,JndGrd,1)
      call LCopy(144,[.false.],0,JfHss,1)
      call ICopy(nSym*16*9,[0],0,JndHss,1)

      ! Overwrite with information in IfGrd, IndGrd, IfHss, and IndHss.

      do iAtom=0,1
        do iCar=0,2
          JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
          do iIrrep=0,nSym-1
            JndGrd(iCar,iAtom,iIrrep) = IndGrd(iCar,iAtom,iIrrep)
          end do
          do jAtom=0,1
            do jCar=0,2
              JfHss(iAtom,iCar,jAtom,jCar) = IfHss(iAtom,iCar,jAtom,jCar)
              do iIrrep=0,nSym-1
                JndHss(iAtom,iCar,jAtom,jCar,iIrrep) = IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
              end do
            end do
          end do
        end do
      end do

      ! Derivatives with respect to the operator is computed via
      ! the translational invariance.

      nnIrrep = nSym
      if (sIrrep) nnIrrep = 1
      do iIrrep=0,nnIrrep-1
        nDisp = IndDsp(kdc+kCnt,iIrrep)
        do iCar=0,2
          iComp = 2**iCar
          if (TF(kdc+kCnt,iIrrep,iComp)) then
            nDisp = nDisp+1

            ! Reset flags for the basis set centers so that we
            ! will explicitly compute the derivatives with
            ! respect to those centers. Activate flag for the
            ! third center so that its derivative will be computed
            ! by the translational invariance.

            JndGrd(iCar,0,iIrrep) = abs(JndGrd(iCar,0,iIrrep))
            JndGrd(iCar,1,iIrrep) = abs(JndGrd(iCar,1,iIrrep))
            JndGrd(iCar,2,iIrrep) = -nDisp
            JfGrd(iCar,0) = .true.
            JfGrd(iCar,1) = .true.
            JfGrd(iCar,2) = .false.
          else
            JndGrd(iCar,2,iIrrep) = 0
          end if
        end do
      end do

      ! The third center is calculated by translational invariance.
      ! This requires the 2nd derivatives on the other centers.

      call LCopy(4,[.false.],0,Tr,1)
      do iCar=0,2
        do jAtom=0,2
          if (jAtom == 2) then
            iStop = iCar
          else
            iStop = 2
          end if
          do jCar=0,iStop
            do iIrrep=0,nSym-1
              if ((JndGrd(iCar,2,iIrrep) /= 0) .and. (JndGrd(jCar,jAtom,iIrrep) /= 0)) then
                JndHss(2,iCar,jAtom,jCar,iIrrep) = -itri(abs(JndGrd(iCar,2,iIrrep)),abs(JndGrd(jCar,jAtom,iIrrep)))

                Tr(2) = .true.
                if (jAtom == 2) then
                  Maxi = max(iCar,jCar)
                  Mini = min(iCar,jCar)
                  JfHss(0,Maxi,0,Mini) = .true.
                  JfHss(1,Maxi,1,Mini) = .true.
                  JfHss(1,iCar,0,jCar) = .true.
                  JfHss(1,jCar,0,iCar) = .true.
                else
                  Maxi = max(iCar,jCar)
                  Mini = min(iCar,jCar)
                  JfHss(jAtom,Maxi,jAtom,Mini) = .true.
                  JfHss(1,iCar,0,jCar) = .true.
                  JfHss(1,jCar,0,iCar) = .true.
                end if
              end if
            end do
          end do
        end do
      end do

      IfG(0) = .true.
      IfG(1) = .true.
      IfG(2) = .false.
      IfG(3) = .false.
      do iCent=0,1
        if (EQ(Coori(1,iCent+1),Coori(1,3))) then
          IfG(iCent) = .false.
          do iCar=0,2
            jfGrd(iCar,iCent) = .false.
            do kCar=0,2
              do kCent=0,3
                jfHss(iCent,iCar,kCent,kCar) = .false.
                jfHss(kCent,kCar,iCent,iCar) = .false.
                do iIrrep=0,nSym-1
                  jndHss(iCent,iCar,kCent,kCar,iIrrep) = 0
                  jndHss(kCent,kCar,iCent,iCar,iIrrep) = 0
                end do
              end do
            end do
            do iIrrep=0,nSym-1
              jndGrd(iCar,iCent,iIrrep) = 0
            end do
          end do
        end if
      end do
      call lCopy(12,[.false.],0,jfgrd,1)

      nFinal = 0
      call Rysg2(iAnga,nRys,nZeta,Array(ipA),Array(ipB),[One],[One],Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,Coori,Coori,CoorAC, &
                 Array(ipArr),nArray,TNAI1,Fake,Cff2D,Array(ipDAO),nDAO,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,mOp,iuvwx,ifg,nFinal, &
                 Indx,.false.,.true.,tr)

    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(final)
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
end if

end subroutine NAHss
