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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine DrvN1_EMB(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: to compute the molecular gradient contribution due to the    *
!         inter-subsystem nuclear repulsion energy                     *
!                                                                      *
! Called from: Alaska                                                  *
!                                                                      *
! Author : F. Aquilante, Geneva, Nov 2010                              *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
integer(kind=iwp) :: iCar, iCnt, iCnttp, iCnttp_B, iComp, iDCRR(0:7), igu, iIrrep, iM1xp, iM2xp, iPrint, iR, iRout, jCnt, jCntMx, &
                     jCnttp, LmbdR, mdc, nCnttp_B, ndc, nDCRR, nDisp
real(kind=wp) :: A(3), B(3), CffM1, CffM2, Cnt0M1, Cnt0M2, Cnt1M1, Cnt1M2, df_dr, dfab, dr_da, fab, Fact, Gam, PreFct, r12, RB(3), &
                 ZA, ZAZB, ZB
logical(kind=iwp) :: EQ, TstFnc
character(len=80) :: Lab
real(kind=wp), allocatable :: Charge_B(:)

if (nIrrep > 1) then
  call WarningMessage(2,'Error in DrvN1_Emb')
  write(u6,*) 'Sorry, OFE gradient code does not understand'
  write(u6,*) 'the use of subsystem symmetry!'
  call Quit_OnUserError()
end if
iRout = 33
iPrint = nPrint(iRout)

iIrrep = 0
Temp(:) = Zero

call NameRun('AUXRFIL') ! switch RUNFILE name

call mma_allocate(Charge_B,nCnttp,label='B-Charges')
call Get_dArray('Nuclear charge',Charge_B,nCnttp)

call NameRun('#Pop')    ! switch back to old RUNFILE name

ZA = One
iCnttp = 1
do while ((iCnttp <= nCnttp) .and. (ZA > Zero))
  ZA = dbsc(iCnttp)%Charge
  iCnttp = iCnttp+1
end do
iCnttp_B = iCnttp-1  ! start of atoms of subsystem B

if (iCnttp_B == 1) then ! subsystem B comes first
  ZB = Zero
  nCnttp_B = 1
  do while ((nCnttp_B <= nCnttp) .and. (ZB == Zero))
    ZB = dbsc(nCnttp_B)%Charge
    nCnttp_B = nCnttp_B+1
  end do
  nCnttp_B = nCnttp_B-1  ! end of atoms of subsystem B
else
  nCnttp_B = nCnttp
end if

mdc = 0
! Loop over centers with the same charge (A-subsystem)
do iCnttp=1,nCnttp
  if (iCnttp > 1) mdc = mdc+dbsc(iCnttp-1)%nCntr
  ZA = dbsc(iCnttp)%Charge
  if ((iCnttp >= iCnttp_B) .and. (iCnttp <= nCnttp_B) .and. (ZA > Zero)) then
    call WarningMessage(2,'Internal error in DrvN1_Emb')
    write(u6,*) ' Subsystems must come one after the other'
    call Abend()
  end if
  if (ZA == Zero) cycle
  ! Loop over all unique centers of this group (A-subsystem)
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=iCnttp_B,nCnttp_B  ! (B-subsystem)

      if (jCnttp > iCnttp_B) ndc = ndc+dbsc(jCnttp-1)%nCntr

      ZB = Charge_B(jCnttp)

      if (ZB == Zero) cycle
      ZAZB = ZA*ZB
      jCntMx = dbsc(jCnttp)%nCntr
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
          if (EQ(A,RB)) cycle
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          fab = One
          dfab = Zero
          if (dbsc(iCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            do iM1xp=1,dbsc(iCnttp)%nM1
              Gam = dbsc(iCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gam*r12**2))
              Cnt1M1 = Cnt1M1+Gam*(CffM1*exp(-Gam*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            do iM2xp=1,dbsc(iCnttp)%nM2
              Gam = dbsc(iCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gam*r12**2))
              Cnt1M2 = Cnt1M2+Gam*(CffM2*exp(-Gam*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
          end if
          if (dbsc(jCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            do iM1xp=1,dbsc(jCnttp)%nM1
              Gam = dbsc(jCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gam*r12**2))
              Cnt1M1 = Cnt1M1+Gam*(CffM1*exp(-Gam*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            do iM2xp=1,dbsc(jCnttp)%nM2
              Gam = dbsc(jCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gam*r12**2))
              Cnt1M2 = Cnt1M2+Gam*(CffM2*exp(-Gam*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
          end if
          df_dr = (dfab*r12-fab)/r12**2

          if (.not. dbsc(iCnttp)%pChrg) then
            nDisp = IndDsp(mdc+iCnt,iIrrep)
            igu = nIrrep/dc(mdc+iCnt)%nStab
            do iCar=0,2
              dr_dA = (A(iCar+1)-RB(iCar+1))/r12
              iComp = 2**iCar
              if (TstFnc(dc(mdc+iCnt)%iCoSet,iIrrep,iComp,dc(mdc+iCnt)%nStab)) then
                nDisp = nDisp+1
                if (Direct(nDisp)) then
                  Temp(nDisp) = Temp(nDisp)+One/real(igu,kind=wp)*PreFct*dr_dA*df_dr
                end if
              end if
            end do
          end if
        end do
      end do
    end do
  end do
end do
if (iPrint >= 15) then
  Lab = ' OFE Nuclear Repulsion Contribution'
  call PrGrad(Lab,Temp,nGrad,ChDisp)
end if

call mma_deallocate(Charge_B)

call DaXpY_(nGrad,One,Temp,1,Grad,1)

return

end subroutine DrvN1_EMB
