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

subroutine M1Grd( &
#                define _CALLING_
#                include "grd_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of the M1 integrals used  *
!         ECP calculations. The operator is the nuclear attraction     *
!         operator times a s-type gaussian function.                   *
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
!                                                                      *
!             Modified to gradients, December '93 (RL).                *
!***********************************************************************

use Basis_Info
use Center_Info

implicit real*8(A-H,O-Z)
external TNAI1, Fake, Cff2D
#include "Molcas.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"
#include "grd_interface.fh"
! Local variables
real*8 C(3), TC(3), CoorAC(3,2), Coori(3,4), Coora(3,4)
integer iDCRT(0:7), iAnga(4), iuvwx(4), lOp(4), JndGrd(3,4)
logical EQ, JfGrad(3,4)
logical, external :: TF
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 193
iPrint = nPrint(iRout)

nRys = nHer

if (iPrint >= 49) then
  call RecPrt(' In M1Grd: A',' ',A,1,3)
  call RecPrt(' In M1Grd: RB',' ',RB,1,3)
  call RecPrt(' In M1Grd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In M1Grd: P',' ',P,nZeta,3)
  write(6,*) ' In M1Grd: la,lb=',' ',la,lb
end if

! Allocate Scratch for primitives and work area for HRR

ip = 1
ipA = ip
ip = ip+nZeta
ipB = ip
ip = ip+nZeta
ipDAO = ip
ip = ip+nZeta*nElem(la)*nElem(lb)
ipK = ip
ip = ip+nZeta
ipZ = ip
ip = ip+nZeta
ipZI = ip
ip = ip+nZeta
ipPx = ip
ip = ip+nZeta
ipPy = ip
ip = ip+nZeta
ipPz = ip
ip = ip+nZeta
if (ip-1 > nArr*nZeta) then
  write(6,*) ' ip-1 > nArr*nZeta (M1 section)'
  write(6,*) ' nArr,nZeta=',nArr,nZeta
  call Abend()
end if
nArray = nArr*nZeta-ip+1

iIrrep = 0
iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
call dcopy_(3,A,1,Coora(1,1),1)
call dcopy_(3,RB,1,Coora(1,2),1)
if (la >= lb) then
  call dcopy_(3,A,1,CoorAC(1,1),1)
else
  call dcopy_(3,RB,1,CoorAC(1,1),1)
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

! Loop over nuclear centers.

kdc = 0
do kCnttp=1,nCnttp
  if (.not. dbsc(kCnttp)%ECP) Go To 111
  if (dbsc(kCnttp)%nM1 == 0) Go To 111
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab

    do lDCRT=0,nDCRT-1
      lOp(3) = NrOpr(iDCRT(lDCRT))
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      ! Branch out if one-center integral
      if (EQ(A,RB) .and. EQ(A,TC)) Go To 102
      if (iPrint >= 99) call RecPrt(' In M1Grd: TC',' ',TC,1,3)
      call dcopy_(3,A,1,Coora(1,1),1)
      call dcopy_(3,RB,1,Coora(1,2),1)
      call dcopy_(6,Coora(1,1),1,Coori(1,1),1)
      if ((.not. EQ(A,RB)) .or. (.not. EQ(A,TC))) then
        Coori(1,1) = Coori(1,1)+One
        !Coora(1,1) = Coora(1,1)+One
      end if
      call dcopy_(3,TC,1,CoorAC(1,2),1)
      call dcopy_(3,TC,1,Coori(1,3),1)
      call dcopy_(3,TC,1,Coori(1,4),1)
      call dcopy_(3,TC,1,Coora(1,3),1)
      call dcopy_(3,TC,1,Coora(1,4),1)

      do iM1xp=1,dbsc(kCnttp)%nM1
        Gamma = dbsc(kCnttp)%M1xp(iM1xp)

        call ICopy(6,IndGrd,1,JndGrd,1)
        do i=1,3
          do j=1,2
            JfGrad(i,j) = IfGrad(i,j)
          end do
        end do

        ! Derivatives with respect to the operator is computed
        ! via the translational invariance.
        ! Some extra care is needed here due to that Rys2Dg will
        ! try to avoid some of the work.

        nDisp = IndDsp(kdc+kCnt,iIrrep)
        do iCar=0,2
          ! No direct assembly of contribution from the operat.
          JfGrad(iCar+1,3) = .false.
          JndGrd(iCar+1,3) = 0
          iCmp = 2**iCar
          if (TF(kdc+kCnt,iIrrep,iCmp) .and. (.not. dbsc(kCnttp)%pChrg)) then
            ! Displacement is symmetric
            nDisp = nDisp+1
            if (Direct(nDisp)) then
              ! Reset flags for the basis set centers so that
              ! we will explicitly compute the derivatives
              ! with respect to those centers. Activate flag
              ! for the third center so that its derivative
              ! will be computed by the translational
              ! invariance.
              JfGrad(iCar+1,1) = .true.
              JfGrad(iCar+1,2) = .true.
              if ((A(iCar+1) /= TC(iCar+1)) .and. (RB(iCar+1) /= TC(iCar+1))) then
                ! Three center case
                JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
                JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
                JndGrd(iCar+1,3) = -nDisp
                JfGrad(iCar+1,1) = .true.
                JfGrad(iCar+1,2) = .true.
              else if ((A(iCar+1) == TC(iCar+1)) .and. (RB(iCar+1) /= TC(iCar+1))) then
                ! Two center case
                JndGrd(iCar+1,1) = -abs(JndGrd(iCar+1,1))
                JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
                JfGrad(iCar+1,1) = .false.
                JfGrad(iCar+1,2) = .true.
              else if ((A(iCar+1) /= TC(iCar+1)) .and. (RB(iCar+1) == TC(iCar+1))) then
                ! Two center case
                JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
                JndGrd(iCar+1,2) = -abs(JndGrd(iCar+1,2))
                JfGrad(iCar+1,1) = .true.
                JfGrad(iCar+1,2) = .false.
              else
                ! One center case
                JndGrd(iCar+1,1) = 0
                JndGrd(iCar+1,2) = 0
                JfGrad(iCar+1,1) = .false.
                JfGrad(iCar+1,2) = .false.
              end if
            end if
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

        ! Modify the original basis. Observe that
        ! simplification due to A=B are not valid for the
        ! exponent index, eq. P-A=/=0.

        do iZeta=1,nZeta
          PTC2 = (P(iZeta,1)-TC(1))**2+(P(iZeta,2)-TC(2))**2+(P(iZeta,3)-TC(3))**2
          Tmp0 = Zeta(iZeta)+Gamma
          Tmp1 = exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
          Array(ipK+iZeta-1) = rKappa(iZeta)*Tmp1
          Array(ipZ+iZeta-1) = Tmp0
          Array(ipZI+iZeta-1) = One/Tmp0
          Array(ipPx+iZeta-1) = (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
          Array(ipPy+iZeta-1) = (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
          Array(ipPz+iZeta-1) = (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
        end do

        ! Modify the density matrix with the prefactor

        Fact = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M1cf(iM1xp)*(dble(nStabM)/dble(LmbdT))*Two*Pi
        nDAO = nElem(la)*nElem(lb)
        do iDAO=1,nDAO
          do iZeta=1,nZeta
            Fac = Fact*Array(ipK+iZeta-1)*Array(ipZI+iZeta-1)
            ipDAOt = nZeta*(iDAO-1)+iZeta-1+ipDAO
            Array(ipDAOt) = Fac*DAO(iZeta,iDAO)
          end do
        end do
        if (iPrint >= 99) then
          write(6,*) ' Charge=',dbsc(kCnttp)%Charge
          write(6,*) ' Fact=',Fact
          write(6,*) ' IndGrd=',IndGrd
          write(6,*) ' JndGrd=',JndGrd
          call RecPrt('DAO*Fact',' ',Array(ipDAO),nZeta,nDAO)
        end if

        ! Compute integrals with the Rys quadrature.

        call Rysg1(iAnga,nRys,nZeta,Array(ipA),Array(ipB),[One],[One],Array(ipZ),Array(ipZI),nZeta,[One],[One],1,Array(ipPx), &
                   nZeta,TC,1,Coori,Coora,CoorAC,Array(ip),nArray,TNAI1,Fake,Cff2D,Array(ipDAO),nDAO,Grad,nGrad,JfGrad,JndGrd,lOp, &
                   iuvwx)

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
  call Unused_real_array(final)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
end if

end subroutine M1Grd
