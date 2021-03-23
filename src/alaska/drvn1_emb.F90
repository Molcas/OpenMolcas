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

use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "SysDef.fh"
#include "print.fh"
#include "real.fh"
#include "disp.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "status.fh"
real*8 A(3), B(3), RB(3), Grad(nGrad), Temp(nGrad)
integer iDCRR(0:7)
logical EQ, TstFnc
character*16 NamRfil
character*80 Lab

!***********************************************************************
!     Statement function for Charges of subsystem B
Charge_B(i) = Work(ip_ChargeB+i-1)
!***********************************************************************

if (nIrrep > 1) then
  call WarningMessage(2,'Error in DrvN1_Emb')
  write(6,*) 'Sorry, OFE gradient code does not understand'
  write(6,*) 'the use of subsystem symmetry!'
  call Quit_OnUserError()
end if
iRout = 33
iPrint = nPrint(iRout)

iIrrep = 0
call dcopy_(nGrad,[Zero],0,Temp,1)

call Get_NameRun(NamRfil) ! save the old RUNFILE name
call NameRun('AUXRFIL')   ! switch RUNFILE name

call GetMem('B-Charges','Allo','Real',ip_ChargeB,nCnttp)
call Get_dArray('Nuclear charge',Work(ip_ChargeB),nCnttp)

call NameRun(NamRfil)   ! switch back to old RUNFILE name

ZA = 1.0d0
iCnttp = 1
do while ((iCnttp <= nCnttp) .and. (ZA > 0.0d0))
  ZA = dbsc(iCnttp)%Charge
  iCnttp = iCnttp+1
end do
iCnttp_B = iCnttp-1  ! start of atoms of subsystem B

if (iCnttp_B == 1) then ! subsystem B comes first
  ZB = 0.0d0
  nCnttp_B = 1
  do while ((nCnttp_B <= nCnttp) .and. (ZB == 0.0d0))
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
  ZA = dbsc(iCnttp)%Charge
  if ((iCnttp >= iCnttp_B) .and. (iCnttp <= nCnttp_B) .and. (ZA > 0.0d0)) then
    call WarningMessage(2,'Internal error in DrvN1_Emb')
    write(6,*) ' Subsystems must come one after the other'
    call Abend()
  end if
  if (ZA == Zero) Go To 101
  ! Loop over all unique centers of this group (A-subsystem)
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=iCnttp_B,nCnttp_B  ! (B-subsystem)

      ZB = Charge_B(jCnttp)

      if (ZB == Zero) Go To 201
      ZAZB = ZA*ZB
      jCntMx = dbsc(jCnttp)%nCntr
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
          if (EQ(A,RB)) Go To 301
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
              Gamma = dbsc(iCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12**2))
              Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            do iM2xp=1,dbsc(iCnttp)%nM2
              Gamma = dbsc(iCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12**2))
              Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12**2))
            end do
            fab = fab+r12*Cnt0M2
            dfab = dfab+(Cnt0M2-Two*r12**2*Cnt1M2)
          end if
          if (dbsc(jCnttp)%ECP) then
            ! Add contribution from M1 operator
            Cnt0M1 = Zero
            Cnt1M1 = Zero
            do iM1xp=1,dbsc(jCnttp)%nM1
              Gamma = dbsc(jCnttp)%M1xp(iM1xp)
              CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
              Cnt0M1 = Cnt0M1+(CffM1*exp(-Gamma*r12**2))
              Cnt1M1 = Cnt1M1+Gamma*(CffM1*exp(-Gamma*r12**2))
            end do
            fab = fab+Cnt0M1
            dfab = dfab-Two*r12*Cnt1M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            do iM2xp=1,dbsc(jCnttp)%nM2
              Gamma = dbsc(jCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gamma*r12**2))
              Cnt1M2 = Cnt1M2+Gamma*(CffM2*exp(-Gamma*r12**2))
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
                  Temp(nDisp) = Temp(nDisp)+One/dble(igu)*PreFct*dr_dA*df_dr
                end if
              end if
            end do
          end if
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
if (iPrint >= 15) then
  Lab = ' OFE Nuclear Repulsion Contribution'
  call PrGrad(Lab,Temp,nGrad,ChDisp,5)
end if

call GetMem('B-Charges','Free','Real',ip_ChargeB,nCnttp)

call DaXpY_(nGrad,One,Temp,1,Grad,1)

return

end subroutine DrvN1_EMB
