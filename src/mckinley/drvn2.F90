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

subroutine DrvN2(Hess,nGrad)
!***********************************************************************
!                                                                      *
! Object: to compute the molecular gradient contribution due to the    *
!         nuclear repulsion energy.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October 1991                                             *
!             Anders Bernhardsson, Dept. of Theoretical Chemistry,     *
!             University of Lund, SWEDEN                               *
!             September 1995                                           *
!***********************************************************************

use Basis_Info
use Center_Info
use PCM_arrays
use Symmetry_Info, only: nIrrep, iChTbl

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "disp.fh"
#include "disp2.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
real*8 A(3), B(3), RB(3), Hess(nGrad*(nGrad+1)/2), prmt(0:7), C(3), D(3), SD(3)
integer iDCRR(0:7), IndGrd(3,2,0:7), ii(2), iStb(0:7), iDCRS(0:7), IndHss(2,3,2,3,0:7), nop(2), kop(2)
logical EQ, NoLoop
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
real*8, allocatable :: Pcmhss(:), Der1(:), DerDM(:), Temp(:)
logical, external :: TF
! Statement Function
xPrmt(i,j) = Prmt(iand(i,j))
iTri(i1,i2) = max(i1,i2)*(max(i1,i2)-1)/2+min(i1,i2)

!                                                                      *
!***********************************************************************
!                                                                      *
!iRout = 33
!iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the nuclear repulsion contributions
!                                                                      *
!***********************************************************************
!                                                                      *
nHess = nGrad*(nGrad+1)/2
call dcopy_(nHess,[Zero],0,Hess,1)

mdc = 0
! Loop over centers with the same change
do iCnttp=1,nCnttp
  ZA = dbsc(iCnttp)%Charge
  if (ZA == Zero) Go To 101
  ! Loop over all unique centers of this group
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=1,iCnttp
      ZB = dbsc(jCnttp)%Charge
      if (ZB == Zero) Go To 201
      ZAZB = ZA*ZB
      jCntMx = dbsc(jCnttp)%nCntr
      if (iCnttp == jCnttp) jCntMx = iCnt
      do jCnt=1,jCntMx
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        Fact = One
        ! Factor due to resticted summation
        if (EQ(A,B)) Fact = Half

        ! Find the DCR for the two centers

        call DCR(LmbdR,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        PreFct = Fact*ZAZB*dble(nIrrep)/dble(LmbdR)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp(1) = NrOpr(0)
          nOp(2) = NrOpr(iDCRR(iR))
          kop(1) = 0
          kop(2) = iDCRR(iR)
          if (EQ(A,RB)) Go To 301
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          fab = One
          dfab = Zero
          ddfab = Zero

          if (dbsc(iCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            Cnt2M1 = Zero
            do iM1xp=1,dbsc(iCnttp)%nM1
              Gamma = dbsc(iCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12**2))
              Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12**2))
              Cnt2M1 = Cnt2M1+Gamma**2*(CffM1*exp(-Gamma*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ddfab = -Two*Cnt1M1+Four*r12**2*Cnt2M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            Cnt2M2 = Zero
            do iM2xp=1,dbsc(iCnttp)%nM2
              Gamma = dbsc(iCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12**2))
              Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12**2))
              Cnt2M2 = Cnt2M2+Gamma**2*(CffM2*exp(-Gamma*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+Cnt0M2-Two*r12**2*Cnt1M2
            ddfab = ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
          end if
          if (dbsc(jCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            Cnt2M1 = Zero
            do iM1xp=1,dbsc(jCnttp)%nM1
              Gamma = dbsc(jCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12**2))
              Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12**2))
              Cnt2M1 = Cnt2M1+Gamma**2*(CffM1*exp(-Gamma*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ddfab = -Two*Cnt1M1+Four*r12**2*Cnt2M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            Cnt2M2 = Zero
            do iM2xp=1,dbsc(jCnttp)%nM2
              Gamma = dbsc(jCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12**2))
              Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12**2))
              Cnt2M2 = Cnt2M2+Gamma**2*(CffM2*exp(-Gamma*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+Cnt0M2-Two*r12**2*Cnt1M2
            ddfab = ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
          end if

          df_dr = (dfab*r12-fab)/r12**2
          d2f_dr2 = ((ddfab*r12)*r12**2-(dfab*r12-fab)*Two*r12)/r12**4

          call ICopy(nirrep*36,[0],0,Indhss,1)
          call ICopy(nirrep*6,[0],0,indgrd,1)

          ! Determine which displacement in all IR's, each center is associated with

          nnIrrep = nIrrep
          if (sIrrep) nnIrrep = 1

          do iIrrep=0,nnIrrep-1
            nDisp1 = IndDsp(mdc+iCnt,iIrrep)
            nDisp2 = IndDsp(ndc+jCnt,iIrrep)
            do iCar=0,2
              iComp = 2**iCar
              if (TF(mdc+iCnt,iIrrep,iComp)) then
                nDisp1 = nDisp1+1
                IndGrd(iCar+1,1,iIrrep) = nDisp1
              else
                IndGrd(iCar+1,1,iIrrep) = 0
              end if
              iComp = 2**iCar
              if (TF(ndc+jCnt,iIrrep,iComp)) then
                nDisp2 = nDisp2+1
                IndGrd(iCar+1,2,iIrrep) = nDisp2
              else
                IndGrd(iCar+1,2,iIrrep) = 0
              end if
            end do ! iCar
          end do   ! iIrrep

          ! Determine index for each 2nd derivative

          do iIrrep=0,nnIrrep-1
            do iAtom=1,2
              do iCar=1,3
                do jAtom=1,iAtom
                  jCar_Max = 3
                  if (iAtom == jAtom) jCar_Max = iCar
                  do jCar=1,jCar_Max
                    if ((IndGrd(iCar,iAtom,iIrrep) > 0) .and. (IndGrd(jCar,jAtom,iIrrep) > 0)) then

                      IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = iTri(IndGrd(iCar,iAtom,iIrrep),IndGrd(jCar,jAtom,iIrrep))

                    else

                      IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0

                    end if
                  end do ! jCar
                end do   ! jAtom
              end do     ! iCar
            end do       ! iAtom
          end do         ! iIrrep

          ii(1) = dc(mdc+icnt)%nStab
          ii(2) = dc(ndc+jcnt)%nStab

          do iIrrep=0,nnIrrep-1
            do iCent=1,2
              do jCent=1,iCent
                do iCar=1,3
                  jCar_Max = 3
                  if (iCent == jCent) jCar_Max = iCar
                  do jCar=1,jCar_Max
                    iCh1 = 2**(iCar-1)
                    iCh2 = 2**(jCar-1)
                    g = dble(iChTbl(iIrrep,nOp(icent)))*xPrmt(kOp(icent),iCh1)*dble(ii(icent))/dble(nIrrep)
                    g = g*dble(iChTbl(iIrrep,nOp(jcent)))*xPrmt(kOp(jcent),iCh2)*dble(ii(jcent))/dble(nIrrep)
                    g = g*(-One)**(icent+jcent)
                    if ((iCent /= jCent) .and. (iCar == jCar) .and. &
                        (abs(indgrd(iCar,iCent,iIrrep)) == abs(indgrd(jCar,jCent,iIrrep)))) then
                      ps = Two
                    else
                      ps = One
                    end if

                    Index = indHss(iCent,iCar,jCent,jCar,iIrrep)
                    if (index /= 0) then
                      dr_dAi = (A(iCar)-RB(iCar))/r12
                      dr_dAj = (A(jCar)-RB(jCar))/r12
                      d2r_dAidAj = -(A(iCar)-RB(iCar))*dr_dAj
                      if (iCar == jCar) d2r_dAidAj = d2r_dAidAj+r12
                      d2r_dAidAj = d2r_dAidAj/r12**2
                      Hess(Index) = Hess(index)+g*PreFct*ps*(d2r_dAidAj*df_dr+dr_dAi*dr_dAj*d2f_dr2)
                    end if
                  end do ! jCar
                end do   ! iCar
              end do     ! jCent
            end do       ! iCent
          end do         ! iIrrep

          !call triprt(' ',' ',Hess,ldisp(0))
301       continue
        end do

      end do
201   continue
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
  end do
101 continue
  mdc = mdc+dbsc(iCnttp)%nCntr
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! PCM contributions
!                                                                      *
!***********************************************************************
!                                                                      *
if (PCM) then

  ! We will have three contributions here.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! 1) Process the contribution
  !
  !    Sum(i) q_i V_in^xy

  ! Loop over tiles

  do iTs=1,nTs
    ZA = PCM_SQ(1,iTs)+PCM_SQ(2,iTs)
    NoLoop = ZA == Zero
    if (NoLoop) Go To 112
    ZA = ZA/dble(nIrrep)
    A(1:3) = PCMTess(1:3,iTs)

    ! Tile only stabilized by the unit operator

    nStb = 1
    iStb(0) = 0

    ndc = 0
    do jCnttp=1,nCnttp
      ZB = dbsc(jCnttp)%Charge
      if (ZB == Zero) Go To 212
      if (dbsc(jCnttp)%pChrg) Go To 212
      if (dbsc(jCnttp)%Frag) Go To 212
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        PreFct = ZAZB*dble(nIrrep)/dble(LmbdR)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp(1) = NrOpr(0)
          nOp(2) = NrOpr(iDCRR(iR))
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          fab = One
          dfab = Zero
          ddfab = Zero
          if (dbsc(jCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            Cnt2M1 = Zero
            do iM1xp=1,dbsc(jCnttp)%nM1
              Gamma = dbsc(jCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12**2))
              Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12**2))
              Cnt2M1 = Cnt2M1+Gamma**2*(CffM1*exp(-Gamma*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ddfab = -Two*Cnt1M1+Four*r12**2*Cnt2M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            Cnt2M2 = Zero
            do iM2xp=1,dbsc(jCnttp)%nM2
              Gamma = dbsc(jCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12**2))
              Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12**2))
              Cnt2M2 = Cnt2M2+Gamma**2*(CffM2*exp(-Gamma*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+Cnt0M2-Two*r12**2*Cnt1M2
            ddfab = ddfab-Six**r12*Cnt1M2+Four*r12*Three*Cnt2M2
          end if

          df_dr = (dfab*r12-fab)/r12**2
          d2f_dr2 = ((ddfab*r12)*r12**2-(dfab*r12-fab)*Two*r12)/r12**4

          call ICopy(nirrep*36,[0],0,Indhss,1)
          call ICopy(nirrep*6,[0],0,indgrd,1)

          ! Determine which displacement in all IRs, each center is associated with

          nnIrrep = nIrrep
          if (sIrrep) nnIrrep = 1

          do iIrrep=0,nnIrrep-1
            nDisp1 = IndDsp(ndc+jCnt,iIrrep)
            do iCar=0,2
              iComp = 2**iCar
              if (TF(ndc+jCnt,iIrrep,iComp)) then
                nDisp1 = nDisp1+1
                IndGrd(iCar+1,1,iIrrep) = nDisp1
              else
                IndGrd(iCar+1,1,iIrrep) = 0
              end if
            end do ! iCar
          end do   ! iIrrep

          ! Determine index for each 2nd derivative

          ! Note that each term is only associated with one basis set center.

          iAtom = 1
          jAtom = 1
          do iIrrep=0,nnIrrep-1
            do iCar=1,3
              jCar_Max = iCar
              do jCar=1,jCar_Max
                if ((IndGrd(iCar,iAtom,iIrrep) > 0) .and. (IndGrd(jCar,jAtom,iIrrep) > 0)) then

                  IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = iTri(IndGrd(iCar,iAtom,iIrrep),IndGrd(jCar,jAtom,iIrrep))

                else

                  IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0

                end if
              end do ! jCar
            end do   ! iCar
          end do     ! iIrrep

          ii(1) = dc(ndc+jcnt)%nStab

          iCent = 1
          jCent = 1
          do iIrrep=0,nnIrrep-1
            do iCar=1,3
              jCar_Max = iCar
              do jCar=1,jCar_Max
                iCh1 = 2**(iCar-1)
                iCh2 = 2**(jCar-1)
                g = dble(iChTbl(iIrrep,nOp(icent)))*xPrmt(kOp(icent),iCh1)*dble(ii(icent))/dble(nIrrep)
                g = g*dble(iChTbl(iIrrep,nOp(jcent)))*xPrmt(kOp(jcent),iCh2)*dble(ii(jcent))/dble(nIrrep)
                g = g*(-One)**(icent+jcent)

                Index = indHss(iCent,iCar,jCent,jCar,iIrrep)
                if (Index /= 0) then
                  dr_dAi = (A(iCar)-RB(iCar))/r12
                  dr_dAj = (A(jCar)-RB(jCar))/r12
                  d2r_dAidAj = -(A(iCar)-RB(iCar))*dr_dAj
                  if (iCar == jCar) d2r_dAidAj = d2r_dAidAj+r12
                  d2r_dAidAj = d2r_dAidAj/r12**2
                  Hess(Index) = Hess(Index)+g*PreFct*(d2r_dAidAj*df_dr+dr_dAi*dr_dAj*d2f_dr2)
                end if
              end do ! jCar
            end do   ! iCar
          end do     ! iIrrep

          !call TriPrt(' ',' ',Hess,ldisp(0))

        end do ! End loop over DCR operators, iR

      end do ! End over centers, jCnt
212   continue
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do ! End over basis set types, jCnttp
112 continue
  end do ! End of tiles
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! 2) Process the contribution
  !
  !    Sum(i,j) V_i,n^x Q_ij V_j,n^y

  do iTs=1,nTs
    A(1:3) = PCMTess(1:3,iTs)

    ! Tile only stabilized by the unit operator

    nStb = 1
    iStb(0) = 0

    do jTs=1,iTs
      Fact = Two
      if (jTs == iTs) Fact = One
      Q_ij = PCMDM(iTs,jTs)
      NoLoop = Q_ij == Zero
      if (NoLoop) Go To 122
      C(1:3) = PCMTess(1:3,jTs)

      ! Loop over the basis functions

      mdc = 0
      do iCnttp=1,nCnttp
        ZA = dbsc(iCnttp)%Charge
        if (ZA == Zero) Go To 222
        if (dbsc(iCnttp)%pChrg) Go To 222
        if (dbsc(iCnttp)%Frag) Go To 222

        do iCnt=1,dbsc(iCnttp)%nCntr
          B(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

          ! Find the DCR for the two centers (

          call DCR(LmbdR,iStb,nStb,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,iDCRR,nDCRR)

          PreFct_AB = dble(nIrrep)/dble(LmbdR)
          do iR=0,nDCRR-1
            call OA(iDCRR(iR),B,RB)
            nOp(1) = NrOpr(iDCRR(iR))
            r12_AB = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)
            fab = One
            dfab = Zero
            if (dbsc(iCnttp)%ECP) then
              ! Add contribution from M1 operator
              Cnt0M1 = Zero
              Cnt1M1 = Zero
              do iM1xp=1,dbsc(iCnttp)%nM1
                Gamma = dbsc(iCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
                Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12_AB**2))
                Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12_AB**2))
              end do
              fab = fab+Cnt0M1
              dfab = dfab-Two*r12_AB*Cnt1M1
              ! Add contribution from M2 operator
              Cnt0M2 = Zero
              Cnt1M2 = Zero
              do iM2xp=1,dbsc(iCnttp)%nM2
                Gamma = dbsc(iCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
                Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12_AB**2))
                Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12_AB**2))
              end do
              fab = fab+r12_AB*Cnt0M2
              dfab = dfab+Cnt0M2-Two*r12_AB**2*Cnt1M2
            end if
            df_dr_AB = (dfab*r12_AB-fab)/r12_AB**2

            ndc = 0
            do jCnttp=1,nCnttp
              ZB = dbsc(jCnttp)%Charge
              if (ZB == Zero) Go To 232
              if (dbsc(jCnttp)%pChrg) Go To 232
              if (dbsc(jCnttp)%Frag) Go To 232

              do jCnt=1,dbsc(jCnttp)%nCntr
                D(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

                ! Find the DCR for the two centers

                call DCR(LmbdS,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRS,nDCRS)

                PreFct_CD = dble(nIrrep)/dble(LmbdS)
                do iS=0,nDCRS-1
                  call OA(iDCRS(iS),D,SD)
                  nOp(2) = NrOpr(iDCRS(iS))
                  r12_CD = sqrt((C(1)-SD(1))**2+(C(2)-SD(2))**2+(C(3)-SD(3))**2)

                  fcd = One
                  dfcd = Zero
                  if (dbsc(jCnttp)%ECP) then
                    ! Add contribution from M1 operator
                    Cnt0M1 = Zero
                    Cnt1M1 = Zero
                    do iM1xp=1,dbsc(jCnttp)%nM1
                      Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                      CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                      Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12_CD**2))
                      Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12_CD**2))
                    end do
                    fcd = fcd+Cnt0M1
                    dfcd = dfcd-Two*r12_CD*Cnt1M1
                    ! Add contribution from M2 operator
                    Cnt0M2 = Zero
                    Cnt1M2 = Zero
                    do iM2xp=1,dbsc(jCnttp)%nM2
                      Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                      CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                      Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12_CD**2))
                      Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12_CD**2))
                    end do
                    fcd = fcd+r12_CD*Cnt0M2
                    dfcd = dfcd+Cnt0M2-Two*r12_CD**2*Cnt1M2
                  end if
                  df_dr_CD = (dfcd*r12_CD-fcd)/r12_CD**2

                  call ICopy(nirrep*36,[0],0,Indhss,1)
                  call ICopy(nirrep*6,[0],0,indgrd,1)

                  ! Determine which displacement in all IR's, each center is associated with

                  nnIrrep = nIrrep
                  if (sIrrep) nnIrrep = 1

                  do iIrrep=0,nnIrrep-1
                    nDisp1 = IndDsp(mdc+iCnt,iIrrep)
                    nDisp2 = IndDsp(ndc+jCnt,iIrrep)
                    do iCar=0,2
                      iComp = 2**iCar

                      if (TF(mdc+iCnt,iIrrep,iComp)) then
                        nDisp1 = nDisp1+1
                        IndGrd(iCar+1,1,iIrrep) = nDisp1
                      else
                        IndGrd(iCar+1,1,iIrrep) = 0
                      end if

                      if (TF(ndc+jCnt,iIrrep,iComp)) then
                        nDisp2 = nDisp2+1
                        IndGrd(iCar+1,2,iIrrep) = nDisp2
                      else
                        IndGrd(iCar+1,2,iIrrep) = 0
                      end if
                    end do ! iCar
                  end do   ! iIrrep

                  ! Determine index for each 2nd derivative

                  do iIrrep=0,nnIrrep-1
                    do iAtom=1,2
                      do iCar=1,3
                        do jAtom=1,iAtom
                          jCar_Max = 3
                          if (iAtom == jAtom) jCar_Max = iCar
                          do jCar=1,jCar_Max
                            if ((IndGrd(iCar,iAtom,iIrrep) > 0) .and. (IndGrd(jCar,jAtom,iIrrep) > 0)) then

                              IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = iTri(IndGrd(iCar,iAtom,iIrrep),IndGrd(jCar,jAtom,iIrrep))

                            else

                              IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0

                            end if
                          end do ! jCar
                        end do   ! jAtom
                      end do     ! iCar
                    end do       ! iAtom
                  end do         ! iIrrep

                  ii(1) = dc(mdc+icnt)%nStab
                  ii(2) = dc(ndc+jcnt)%nStab

                  ! Note that we have two different cases here, depending on if
                  ! iTs=jTs or not!
                  ! For iTs=jTs and iCent == jCent we do
                  ! dV_i/dx*dV_i/dy only and exclude dV_i/dy*dV_i/dx since they
                  ! are the same. For iTs /= jTs we need both dV_i/dx*dV_j/dy and
                  ! dV_i/dy*dV_j/dx.

                  do iIrrep=0,nnIrrep-1
                    do iCent=1,2
                      do jCent=1,iCent
                        do iCar=1,3
                          jCar_Max = 3
                          if ((iCent == jCent) .and. (iTs == jTs)) jCar_Max = iCar
                          do jCar=1,jCar_Max
                            iCh1 = 2**(iCar-1)
                            iCh2 = 2**(jCar-1)
                            g = dble(iChTbl(iIrrep,nOp(icent)))*xPrmt(kOp(icent),iCh1)*dble(ii(icent))/dble(nIrrep)
                            g = g*dble(iChTbl(iIrrep,nOp(jcent)))*xPrmt(kOp(jcent),iCh2)*dble(ii(jcent))/dble(nIrrep)
                            g = g*(-One)**(icent+jcent)

                            if ((iCent /= jCent) .and. (iCar == jCar) .and. &
                                (abs(IndGrd(iCar,iCent,iIrrep)) == abs(IndGrd(jCar,jCent,iIrrep)))) then
                              ps = Two
                            else
                              ps = One
                            end if

                            Index = IndHss(iCent,iCar,jCent,jCar,iIrrep)
                            if (Index /= 0) then
                              dr_dB = -(A(iCar)-RB(iCar))/r12_AB
                              dr_dD = -(C(jCar)-SD(jCar))/r12_CD
                              Hess(Index) = Hess(Index)+Fact*g*ps*ZA*ZA*Q_ij*PreFct_AB*dr_dB*df_dr_AB*PreFct_CD*dr_dD*df_dr_CD
                            end if
                          end do ! jCar
                        end do   ! iCar
                      end do     ! jCent
                    end do       ! iCent
                  end do         ! iIrrep

                  !call TriPrt(' ',' ',Hess,ldisp(0))

                end do ! iS

              end do   ! jCnt
232           continue
              ndc = ndc+dbsc(jCnttp)%nCntr
            end do     ! jCnttp

          end do       ! iR

        end do         ! iCnt
222     continue
        mdc = mdc+dbsc(iCnttp)%nCntr
      end do           ! jCnttp
122   continue
    end do             ! jTs
  end do               ! iTs
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Add additional contributions

  nPCMHss = nGrad*nGrad
  call Get_nAtoms_All(nAtoms)
  call mma_allocate(pcmhss,nPCMHss,Label='pcmhss')
  call mma_allocate(Der1,nTs,Label='Der1')
  call mma_allocate(DerDM,nTs*nTs,Label='DerDM')
  call mma_allocate(Temp,nTs*nTs,Label='Temp')
  call Cav_Hss(nAtoms,nGrad,nTs,nS,Eps,PCMSph,PCMiSph,PCM_N,PCMTess,PCM_SQ,PCMDM,Der1,DerDM,Temp,dTes,DPnt,dRad,dCntr,pcmhss)
  call mma_deallocate(pcmhss)
  call mma_deallocate(Der1)
  call mma_deallocate(DerDM)
  call mma_deallocate(Temp)

end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine DrvN2
