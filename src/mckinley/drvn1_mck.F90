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

subroutine DrvN1_mck(Grad,nGrad)
!***********************************************************************
!                                                                      *
! Object: to compute the molecular gradient contribution due to the    *
!         nuclear repulsion energy.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October 1991                                             *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iChBas
use Disp, only: IndDsp
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
integer(kind=iwp) :: iCar, iCnt, iCnttp, iComp, iDCRR(0:7), igu, igv, iIrrep, iR, jCnt, jCntMx, jCnttp, LmbdR, mdc, ndc, nDCRR, &
                     nDisp, nOp
real(kind=wp) :: A(3), B(3), Fact, PreFct, ps, r12, RB(3), ZA, ZAZB
integer(kind=iwp), external :: iPrmt, NrOpr
logical(kind=iwp), external :: EQ, TstFnc

iIrrep = 0
mdc = 0
! Loop over centers with the same change
do iCnttp=1,nCnttp
  if (iCnttp > 1) mdc = mdc+dbsc(iCnttp-1)%nCntr
  ZA = dbsc(iCnttp)%Charge
  if (ZA == Zero) cycle
  ! Loop over all unique centers of this group
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=1,iCnttp
      if (jCnttp > 1) ndc = ndc+dbsc(jCnttp-1)%nCntr
      ZAZB = ZA*dbsc(jCnttp)%Charge
      if (ZAZB == Zero) cycle
      jCntMx = dbsc(jCnttp)%nCntr
      if (iCnttp == jCnttp) jCntMx = iCnt
      do jCnt=1,jCntMx
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        Fact = One
        ! Factor due to resticted summation
        if (EQ(A,B)) Fact = Half

        ! Find the DCR for the two centers

        call DCR(LmbdR,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        PreFct = Fact*ZAZB*real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp = NrOpr(iDCRR(iR))
          if (EQ(A,RB)) cycle
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          nDisp = IndDsp(mdc+iCnt,iIrrep)
          igu = nIrrep/dc(mdc+iCnt)%nStab
          do iCar=0,2
            iComp = 2**iCar
            if (TstFnc(dc(mdc+iCnt)%iCoSet,iIrrep,iComp,dc(mdc+iCnt)%nStab)) then
              nDisp = nDisp+1
              Grad(nDisp) = Grad(nDisp)-One/real(igu,kind=wp)*PreFct*(A(iCar+1)-RB(iCar+1))/(r12**3)
            end if
          end do

          nDisp = IndDsp(ndc+jCnt,iIrrep)
          igv = nIrrep/dc(ndc+jCnt)%nStab
          do iCar=0,2
            iComp = 2**iCar
            if (TstFnc(dc(ndc+jCnt)%iCoSet,iIrrep,iComp,dc(ndc+jCnt)%nStab)) then
              nDisp = nDisp+1
              ps = real(iPrmt(nOp,iChBas(2+iCar)),kind=wp)
              Grad(nDisp) = Grad(nDisp)+ps*One/real(igv,kind=wp)*PreFct*(A(iCar+1)-RB(iCar+1))/(r12**3)
            end if
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine DrvN1_mck
