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
! Copyright (C) 1993, Roland Lindh                                     *
!               1993, Per Boussard                                     *
!***********************************************************************

subroutine M2Grd( &
#                define _CALLING_
#                include "grd_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of M2 integrals used in   *
!         ECP calculations. The operator is a s-type gaussian          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Index_Functions, only: nTri_Elem1
use Disp, only: Dirct, IndDsp
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, iAlpha, iBeta, iCar, iCmp, iDCRT(0:7), iIrrep, iM2xp, ipA, ipAxyz, ipB, ipBxyz, ipK, ipPx, ipPy, ipPz, &
                     ipQxyz, iPrint, ipRxyz, ipZ, iRout, iStrt, iuvwx(4), iZeta, j, JndGrd(3,4), kCnt, kCnttp, kdc, lDCRT, LmbdT, &
                     lOp(4), mGrad, mVec, nDAO, nDCRT, nDisp, nip
real(kind=wp) :: C(3), Fact, Factor, Gmma, PTC2, TC(3), Tmp0, Tmp1
logical(kind=iwp) :: ABeq(3), EQ, JfGrad(3,4)
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: TF
#include "print.fh"

#include "macros.fh"
unused_var(ZInv)

iRout = 122
iPrint = nPrint(iRout)

iIrrep = 0
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)
nDAO = nTri_Elem1(la)*nTri_Elem1(lb)

nip = 1
ipA = nip
nip = nip+nZeta
ipB = nip
nip = nip+nZeta
ipAxyz = nip
nip = nip+nZeta*3*nHer*(la+2)
ipBxyz = nip
nip = nip+nZeta*3*nHer*(lb+2)
ipRxyz = nip
nip = nip+nZeta*3*nHer
ipQxyz = nip
nip = nip+nZeta*3*(la+2)*(lb+2)
ipK = nip
nip = nip+nZeta
ipZ = nip
nip = nip+nZeta
ipPx = nip
nip = nip+nZeta
ipPy = nip
nip = nip+nZeta
ipPz = nip
nip = nip+nZeta
if (nip-1 > nArr*nZeta) then
  write(u6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  write(u6,*) ' Abend in M2Grd'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In M2Grd: A',' ',A,1,3)
  call RecPrt(' In M2Grd: RB',' ',RB,1,3)
  call RecPrt(' In M2Grd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M2Grd: Kappa',' ',rKappa,nAlpha,nBeta)
  call RecPrt(' In M2Grd: Zeta',' ',Zeta,nAlpha,nBeta)
  call RecPrt(' In M2Grd: P',' ',P,nZeta,3)
  write(u6,*) ' In M2Grd: la,lb,nHer=',la,lb,nHer
end if

iStrt = ipA
do iBeta=1,nBeta
  call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
  iStrt = iStrt+nAlpha
end do

iStrt = ipB
do iAlpha=1,nAlpha
  call dcopy_(nBeta,Beta,1,Array(iStrt),nAlpha)
  iStrt = iStrt+1
end do

! Loop over nuclear centers

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nM2 == 0) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)
    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    do lDCRT=0,nDCRT-1
      lOp(3) = NrOpr(iDCRT(lDCRT))
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      if (EQ(A,RB) .and. EQ(A,TC)) cycle

      do iM2xp=1,dbsc(kCnttp)%nM2
        Gmma = dbsc(kCnttp)%M2xp(iM2xp)
        if (iPrint >= 99) write(u6,*) ' Gmma=',Gmma

        JndGrd(:,1:2) = IndGrd(:,:)
        do i=1,3
          do j=1,2
            JfGrad(i,j) = IfGrad(i,j)
          end do
        end do

        ! Derivatives with respect to the operator is computed
        ! via the translational invariance.

        nDisp = IndDsp(kdc+kCnt,iIrrep)
        do iCar=0,2
          JfGrad(iCar+1,3) = .false.
          iCmp = 2**iCar
          if (TF(kdc+kCnt,iIrrep,iCmp) .and. (.not. dbsc(kCnttp)%pChrg)) then
            nDisp = nDisp+1
            if (Dirct(nDisp)) then
              ! Reset flags for the basis set centers so that
              ! we will explicitly compute the derivatives
              ! with respect to those centers. Activate flag
              ! for the third center so that its derivative
              ! will be computed by the translational
              ! invariance.
              JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
              JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
              JndGrd(iCar+1,3) = -nDisp
              JfGrad(iCar+1,1) = .true.
              JfGrad(iCar+1,2) = .true.
            else
              JndGrd(iCar+1,3) = 0
            end if
          else
            JndGrd(iCar+1,3) = 0
          end if
        end do
        ! No derivatives with respect to the fourth center.
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

        ! Modify the original basis.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gmma
          Tmp1 = exp(-Zeta(iZeta)*Gmma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gmma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gmma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gmma*TC(3))/Tmp0
        end do
        if (iPrint >= 99) then
          write(u6,*) ' The modified basis set'
          call RecPrt(' In M2Grd: Kappa',' ',Array(ipK),nAlpha,nBeta)
          call RecPrt(' In M2Grd: Zeta',' ',Array(ipZ),nAlpha,nBeta)
          call RecPrt(' In M2Grd: P',' ',Array(ipPx),nZeta,3)
          call RecPrt(' In M2Grd: TC',' ',TC,1,3)
        end if

        ! Compute the cartesian values of the basis functions
        ! angular part

        ABeq(1) = (A(1) == RB(1)) .and. (A(1) == TC(1))
        ABeq(2) = (A(2) == RB(2)) .and. (A(2) == TC(2))
        ABeq(3) = (A(3) == RB(3)) .and. (A(3) == TC(3))
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),nHer,ABeq)
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the contribution from the multipole moment operator

        ABeq(1) = .false.
        ABeq(2) = .false.
        ABeq(3) = .false.
        call CrtCmp(Array(ipZ),Array(ipPx),nZeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),nHer,ABeq)

        ! Compute the cartesian components for the multipole
        ! moment integrals. The integrals are factorized into
        ! components.

        call Assmbl(Array(ipQxyz),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,nZeta,HerW(iHerW(nHer)),nHer)

        ! Combine the cartesian components to the full one
        ! electron integral gradient.

        Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M2cf(iM2xp)*Fact
        call CmbnM2(Array(ipQxyz),nZeta,la,lb,Array(ipZ),Array(ipK),rFinal,Array(ipA),Array(ipB),JfGrad,Factor,mVec)
        if (iPrint >= 99) call RecPrt(' rFinal in M2Grd',' ',rFinal,nZeta*nTri_Elem1(la)*nTri_Elem1(lb),mVec)

        ! Distribute the gradient contributions

        call DistG1X(rFinal,DAO,nZeta,nDAO,mVec,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

      end do

    end do
  end do

end do

return

end subroutine M2Grd
