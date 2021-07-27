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
!      Alpha : exponents of bra gaussians                              *
!      nAlpha: number of primitives (exponents) of bra gaussians       *
!      Beta  : as Alpha but for ket gaussians                          *
!      nBeta : as nAlpha but for the ket gaussians                     *
!      Zeta  : sum of exponents (nAlpha x nBeta)                       *
!      ZInv  : inverse of Zeta                                         *
!      rKappa: gaussian prefactor for the products of bra and ket      *
!              gaussians.                                              *
!      P     : center of new gaussian from the products of bra and ket *
!              gaussians.                                              *
!      Final : array for computed integrals                            *
!      nZeta : nAlpha x nBeta                                          *
!      nComp : number of components in the operator (e.g. dipolmoment  *
!              operator has three components)                          *
!      la    : total angular momentum of bra gaussian                  *
!      lb    : total angular momentum of ket gaussian                  *
!      A     : center of bra gaussian                                  *
!      B     : center of ket gaussian                                  *
!      nRys  : order of Rys- or Hermite-Gauss polynomial               *
!      Array : Auxiliary memory as requested by ECPMem                 *
!      nArr  : length of Array                                         *
!      Ccoor : coordinates of the operator, zero for symmetric oper.   *
!      NOrdOp: Order of the operator                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Basis_Info
use Center_Info
use Her_RW

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"
#include "grd_interface.fh"
! Local variables
real*8 TC(3), C(3)
integer iDCRT(0:7), iuvwx(4), lOp(4), JndGrd(3,4)
logical ABeq(3), JfGrad(3,4), EQ
logical, external :: TF
! Statement function for Cartesian index
nElem(k) = (k+1)*(k+2)/2

iRout = 122
iPrint = nPrint(iRout)

iIrrep = 0
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = kOp(1)
lOp(2) = kOp(2)
nDAO = nElem(la)*nElem(lb)

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
  write(6,*) ' nArr is Wrong! ',nip-1,' > ',nArr*nZeta
  call ErrTra()
  write(6,*) ' Abend in M2Grd'
  call Abend()
end if

if (iPrint >= 49) then
  call RecPrt(' In M2Grd: A',' ',A,1,3)
  call RecPrt(' In M2Grd: RB',' ',RB,1,3)
  call RecPrt(' In M2Grd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M2Grd: Kappa',' ',rKappa,nAlpha,nBeta)
  call RecPrt(' In M2Grd: Zeta',' ',Zeta,nAlpha,nBeta)
  call RecPrt(' In M2Grd: P',' ',P,nZeta,3)
  write(6,*) ' In M2Grd: la,lb,nHer=',la,lb,nHer
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
  if (.not. dbsc(kCnttp)%ECP) Go To 111
  if (dbsc(kCnttp)%nM2 == 0) Go To 111

  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)
    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    do lDCRT=0,nDCRT-1
      lOp(3) = NrOpr(iDCRT(lDCRT))
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      if (EQ(A,RB) .and. EQ(A,TC)) Go To 102

      do iM2xp=1,dbsc(kCnttp)%nM2
        Gamma = dbsc(kCnttp)%M2xp(iM2xp)
        if (iPrint >= 99) write(6,*) ' Gamma=',Gamma

        call ICopy(6,IndGrd,1,JndGrd,1)
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
            if (Direct(nDisp)) then
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
        call ICopy(3,[0],0,JndGrd(1,4),1)
        JfGrad(1,4) = .false.
        JfGrad(2,4) = .false.
        JfGrad(3,4) = .false.
        mGrad = 0
        do iCar=1,3
          do i=1,2
            if (JfGrad(iCar,i)) mGrad = mGrad+1
          end do
        end do
        if (iPrint >= 99) write(6,*) ' mGrad=',mGrad
        if (mGrad == 0) Go To 1011

        ! Modify the original basis.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gamma
          Tmp1 = exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
        end do
        if (iPrint >= 99) then
          write(6,*) ' The modified basis set'
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
        call CmbnM2(Array(ipQxyz),nZeta,la,lb,Array(ipZ),Array(ipK),Final,Array(ipA),Array(ipB),JfGrad,Factor,mVec)
        if (iPrint >= 99) call RecPrt(' Final in M2Grd',' ',Final,nZeta*nElem(la)*nElem(lb),mVec)

        ! Distribute the gradient contributions

        call DistG1X(Final,DAO,nZeta,nDAO,mVec,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

1011    continue
      end do

102   continue
    end do
  end do
111 continue
  kdc = kdc+dbsc(kCnttp)%nCntr

end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(ZInv)
  call Unused_integer_array(lOper)
end if

end subroutine M2Grd
