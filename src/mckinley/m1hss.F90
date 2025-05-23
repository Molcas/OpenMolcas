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

subroutine M1Hss( &
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

use McKinley_global, only: sIrrep
use Index_Functions, only: iTri, nTri_Elem1
use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Disp, only: IndDsp
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "hss_interface.fh"
integer(kind=iwp) :: iAnga(4), iBeta, iCar, iCent, iComp, iDCRT(0:7), iIrrep, ipA, ipAOff, ipArr, ipB, ipBOff, iStop, iuvwx(4), &
                     jAtom, jCar, JndGrd(0:2,0:3,0:7), JndHss(0:3,0:2,0:3,0:2,0:7), kCnt, kCnttp, kdc, lDCRT, LmbdT, Maxi, Mini, &
                     mOp(4), nArray, nDAO, nDCRT, nDisp, nip, nnIrrep, nRys
real(kind=wp) :: C(3), CoorAC(3,2), Coori(3,4), Fact, TC(3)
logical(kind=iwp) :: IfG(0:3), JfGrd(0:2,0:3), JfHss(0:3,0:2,0:3,0:2), Tr(0:3)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ, TF

#include "macros.fh"
unused_var(ZInv)
unused_var(Ccoor)
unused_var(nOrdOp)

!if (iPrint >= 99) then
!  write(u6,*) ' In M1Hss: nArr=',nArr
!end if

nRys = nHer

nip = 1
ipA = nip
nip = nip+nAlpha*nBeta
ipB = nip
nip = nip+nAlpha*nBeta
ipArr = nip
nArray = nArr-nip+1

iIrrep = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
Coori(:,1) = A
Coori(:,2) = RB
if (la >= lb) then
  CoorAC(:,1) = A
else
  CoorAC(:,1) = RB
end if
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
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

! Modify the density matrix with the prefactor

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
!do iDAO=1,nDAO
!  DAO(:,iDAO) = Two*rKappa(:)*Pi*ZInv(:)*DAO(:,iDAO)
!end do
!if (iPrint >= 99) call RecPrt('DAO',' ',DAO,nZeta,nDAO)

! Here we go

! Loop over nuclear centers

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (dbsc(kCnttp)%Charge == Zero) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    do lDCRT=0,nDCRT-1
      JndHss(:,:,:,:,0:nSym-1) = 0
      JndGrd(:,:,0:nSym-1) = 0
      JfHss(:,:,:,:) = .false.
      JfGrd(:,:) = .false.
      Tr(:) = .false.
      mOp(3) = NrOpr(iDCRT(lDCRT))
      mOp(4) = mOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      CoorAC(:,2) = TC
      Coori(:,3) = TC
      Coori(:,4) = TC
      if (EQ(A,TC) .and. EQ(A,RB)) cycle

      ! COPY CNTLR MATRICES

      JfGrd(0:2,0:1) = Ifgrd(0:2,0:1)
      JndGrd(0:2,0:1,0:nSym-1) = IndGrd(0:2,0:1,0:nSym-1)
      JfHss(0:1,0:2,0:1,0:2) = IfHss(0:1,0:2,0:1,0:2)
      JndHss(0:1,0:2,0:1,0:2,0:nSym-1) = IndHss(0:1,0:2,0:1,0:2,0:nSym-1)

      Fact = -dbsc(kCnttp)%Charge*real(nStabM,kind=wp)/real(LmbdT,kind=wp)
      !Array(ipDAO:ipDAO+nZeta*nDAO-1) = Fact*reshape(DAO(:,1:nDAO),[nZeta*nDAO])
      iuvwx(3) = dc(kdc+kCnt)%nStab
      iuvwx(4) = dc(kdc+kCnt)%nStab

      ! Derivatives with respect to the operator are computed via the
      ! translational invariance.

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

            JndGrd(iCar,0:1,iIrrep) = abs(JndGrd(iCar,0:1,iIrrep))
            JndGrd(iCar,2,iIrrep) = -nDisp
            JfGrd(iCar,0:1) = .true.
            JfGrd(iCar,2) = .false.
          else
            JndGrd(iCar,2,iIrrep) = 0
          end if
        end do
      end do

      ! The third center is claculated by translation invarians
      ! This requires the 2nd derivatives on the other centers.

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
                JndHss(2,iCar,jAtom,jCar,iIrrep) = -iTri(abs(JndGrd(iCar,2,iIrrep)),abs(JndGrd(jCar,jAtom,iIrrep)))

                Tr(2) = .true.
                if (jAtom == 2) then
                  Maxi = max(iCar,jCar)
                  Mini = min(iCar,jCar)
                  jfHss(0,Maxi,0,Mini) = .true.
                  jfHss(1,Maxi,1,Mini) = .true.
                  jfHss(1,iCar,0,jCar) = .true.
                  jfHss(1,jCar,0,iCar) = .true.
                else
                  Maxi = max(iCar,jCar)
                  Mini = min(iCar,jCar)
                  jfHss(jAtom,Maxi,jAtom,Mini) = .true.
                  jfHss(1,iCar,0,jCar) = .true.
                  jfHss(1,jCar,0,iCar) = .true.
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
          JfGrd(:,iCent) = .false.
          JfHss(iCent,:,:,:) = .false.
          JfHss(:,:,iCent,:) = .false.
          JndHss(iCent,:,:,:,0:nSym-1) = 0
          JndHss(:,:,iCent,:,0:nSym-1) = 0
          JndGrd(:,iCent,0:nSym-1) = 0
        end if
      end do
      JfGrd(:,:) = .false.

      call M1Kernel(rFinal,Hess,nHess,DAO,nDAO,iAnga,nRys,nZeta,Array(ipA),Array(ipB),Zeta,rKappa,P,TC,Coori,CoorAC,Array(ipArr), &
                    nArray,jfgrd,jndgrd,jfhss,jndhss,ifg,tr,mop,iuvwx,kCnttp,Fact,loper(1),0)

    end do
  end do
end do

return

end subroutine M1Hss
