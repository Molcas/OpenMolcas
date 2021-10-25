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
! Copyright (C) 1991,1995, Roland Lindh                                *
!***********************************************************************

subroutine DrvN1(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: to compute the molecular gradient contribution due to the    *
!         nuclear repulsion energy.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!                                                                      *
!             Modified for ECP's and external electric fields, May '95 *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use PCM_arrays, only: PCM_SQ, PCMTess, MM
use External_Centers, only: iXPolType, nOrd_XF, nXF, XF
use Symmetry_Info, only: iChBas, nIrrep
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"
integer(kind=iwp) :: iCar, iChxyz, iCnt, iCnttp, iComp, iDCRR(0:7), iDum, iFd, igu, igv, iIrrep, iM1xp, iM2xp, ip, iPrint, iR, &
                     iRout, iStb(0:7), iTs, ix, iy, iz, jCnt, jCntMx, jCnttp, jCoSet(8,8), LmbdR, mdc, nCav, ndc, nDCRR, nDisp, &
                     nOp, nStb
real(kind=wp) :: A(3), B(3), CCoMx, CCoMxd, CCoMy, CCoMyd, CCoMz, CCoMzd, CffM1, CffM2, Cnt0M1, Cnt0M2, Cnt1M1, Cnt1M2, DA(3), &
                 DARB, df_dr, dfab, dr_dA, dr_dB, fab, fab0, fab1, fab2, Fact, Gam, PreFct, ps, r12, RB(3), Tempd(3), ZA, ZAZB, ZB
logical(kind=iwp) :: EQ, TstFnc, NoLoop
character(len=80) :: Lab
integer(kind=iwp), external :: iChAtm, iPrmt, NrOpr

iRout = 33
iPrint = nPrint(iRout)
!iPrint = 15

iIrrep = 0

!***********************************************************************
!                                                                      *
!            Compute the nuclear repulsion contribution                *
!                                                                      *
!***********************************************************************

Temp(:) = Zero

mdc = 0
! Loop over centers with the same charge
do iCnttp=1,nCnttp
  if (iCnttp > 1) mdc = mdc+dbsc(iCnttp-1)%nCntr
  if (dbsc(iCnttp)%Frag) then
    ZA = dbsc(iCnttp)%FragCharge
  else
    ZA = dbsc(iCnttp)%Charge
  end if
  if (ZA == Zero) cycle
  ! Loop over all unique centers of this group
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=1,iCnttp
      if (jCnttp > 1) ndc = ndc+dbsc(jCnttp-1)%nCntr
      if (dbsc(jCnttp)%Frag) then
        ZB = dbsc(jCnttp)%FragCharge
      else
        ZB = dbsc(jCnttp)%Charge
      end if
      if (ZB == Zero) cycle
      if (dbsc(iCnttp)%pChrg .and. dbsc(jCnttp)%pChrg) cycle
      if (dbsc(iCnttp)%Frag .and. dbsc(jCnttp)%Frag) cycle
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

        PreFct = Fact*ZAZB*real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp = NrOpr(iDCRR(iR))
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

          if (.not. dbsc(jCnttp)%pChrg) then
            nDisp = IndDsp(ndc+jCnt,iIrrep)
            igv = nIrrep/dc(ndc+jCnt)%nStab
            do iCar=0,2
              dr_dB = -(A(iCar+1)-RB(iCar+1))/r12
              iComp = 2**iCar
              if (TstFnc(dc(ndc+jCnt)%iCoSet,iIrrep,iComp,dc(ndc+jCnt)%nStab)) then
                nDisp = nDisp+1
                if (Direct(nDisp)) then
                  ps = real(iPrmt(nOp,iChBas(2+iCar)),kind=wp)
                  Temp(nDisp) = Temp(nDisp)+ps*One/real(igv,kind=wp)*PreFct*dr_dB*df_dr
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
  Lab = ' The Nuclear Repulsion Contribution'
  call PrGrad(Lab,Temp,nGrad,ChDisp)
end if

call DaXpY_(nGrad,One,Temp,1,Grad,1)

!***********************************************************************
!                                                                      *
!           Compute contribution due to the external field             *
!                                                                      *
!***********************************************************************

if (allocated(XF)) then

  if ((nOrd_XF > 1) .or. (iXPolType > 0)) then
    call WarningMessage(2,'Error in DrvN1')
    write(u6,*) 'Sorry, gradients are not implemented for'
    write(u6,*) 'higher XF than dipoles or for polarisabilities'
    call Quit_OnUserError()
  end if

  Temp(:) = Zero

  iDum = 0
  do iFd=1,nXF
    ZA = XF(4,iFd)
    if (nOrd_XF == 0) then
      DA(1:3) = Zero
    else
      DA(1:3) = XF(5:7,iFd)
    end if
    NoLoop = (ZA == Zero) .and. (DA(1) == Zero) .and. (DA(2) == Zero) .and. (DA(3) == Zero)
    if (NoLoop) cycle
    A(1:3) = XF(1:3,iFd)
    iChxyz = iChAtm(A)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    ndc = 0
    do jCnttp=1,nCnttp
      if (jCnttp > 1) ndc = ndc+dbsc(jCnttp-1)%nCntr
      ZB = dbsc(jCnttp)%Charge
      if (ZB == Zero) cycle
      if (dbsc(jCnttp)%pChrg) cycle
      if (dbsc(jCnttp)%Frag) cycle
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        PreFct = real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp = NrOpr(iDCRR(iR))
          if (EQ(A,RB)) cycle
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)
          DARB = DA(1)*(A(1)-RB(1))+DA(2)*(A(2)-RB(2))+DA(3)*(A(3)-RB(3))

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          fab0 = One
          fab1 = One
          fab2 = Three
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
            fab0 = fab0+Cnt0M1-Two*r12**2*Cnt1M1
            fab1 = fab1+Cnt0M1
            fab2 = fab2+Three*Cnt0M1+Two*r12**2*Cnt1M1
            ! Add contribution from M2 operator
            Cnt0M2 = Zero
            Cnt1M2 = Zero
            do iM2xp=1,dbsc(jCnttp)%nM2
              Gam = dbsc(jCnttp)%M2xp(iM2xp)
              CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
              Cnt0M2 = Cnt0M2+(CffM2*exp(-Gam*r12**2))
              Cnt1M2 = Cnt1M2+Gam*(CffM2*exp(-Gam*r12**2))
            end do
            fab0 = fab0+Two*(r12*Cnt0M2-r12**3*Cnt1M2)
            fab1 = fab1+r12*Cnt0M2
            fab2 = fab2-Two*(r12*Cnt0M2-r12**3*Cnt1M2)
          end if

          nDisp = IndDsp(ndc+jCnt,iIrrep)
          igv = nIrrep/dc(ndc+jCnt)%nStab
          do iCar=0,2
            iComp = 2**iCar
            if (TstFnc(dc(ndc+jCnt)%iCoSet,iIrrep,iComp,dc(ndc+jCnt)%nStab)) then
              nDisp = nDisp+1
              if (Direct(nDisp)) then
                ps = real(iPrmt(nOp,iChBas(2+iCar)),kind=wp)
                Temp(nDisp) = Temp(nDisp)+ps*One/real(igv,kind=wp)*PreFct* &
                              (ZAZB*fab0*(A(iCar+1)-RB(iCar+1))/(r12**3)+ &
                               ZB*(fab1*DA(iCar+1)/(r12**3)-DARB*fab2*(A(iCar+1)-RB(iCar+1))/(r12**5)))
              end if
            end if
          end do   ! End loop over cartesian components, iCar

        end do     ! End loop over DCR operators, iR

      end do       ! End over centers, jCnt
    end do         ! End over basis set types, jCnttp
  end do           ! End of centers of the external field, iFD
  if (iPrint >= 15) then
    Lab = ' The Nuclear External Electric Field Contribution'
    call PrGrad(Lab,Temp,nGrad,ChDisp)
  end if

  call DaXpY_(nGrad,One,Temp,1,Grad,1)
end if

!***********************************************************************
!                                                                      *
!          Compute contributions due to the reaction field             *
!                 KirkWood Model                                       *
!                                                                      *
!***********************************************************************

if (lRF .and. (.not. lLangevin) .and. (.not. PCM)) then
  nCav = (lMax+1)*(lMax+2)*(lMax+3)/6

  ! Get the multipole moments

  call Get_dArray('RCTFLD',MM,nCav*2)
  if (iPrint >= 99) call RecPrt('Total Multipole Moments',' ',MM(1,1),1,nCav)
  if (iPrint >= 99) call RecPrt('Total Electric Field',' ',MM(1,2),1,nCav)

  Temp(:) = Zero

  ip = 0
  do ir=0,lMax
    do ix=ir,0,-1
      do iy=ir-ix,0,-1
        iz = ir-ix-iy
        ip = ip+1
        if (iPrint >= 99) write(u6,*) ' ix,iy,iz=',ix,iy,iz

        mdc = 0
        do iCnttp=1,nCnttp
          if (iCnttp > 1) mdc = mdc+dbsc(iCnttp-1)%nCntr
          if (dbsc(iCnttp)%Charge == Zero) cycle
          if (dbsc(iCnttp)%Frag) cycle
          ZA = dbsc(iCnttp)%Charge
          if (iPrint >= 99) then
            write(u6,*) ' Charge=',ZA
            call RecPrt(' Centers',' ',dbsc(iCnttp)%Coor,3,dbsc(iCnttp)%nCntr)
          end if
          do iCnt=1,dbsc(iCnttp)%nCntr
            A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

            if (ix == 0) then
              CCoMx = One
              CCoMxd = Zero
            else if (ix == 1) then
              CCoMx = A(1)
              CCoMxd = One
            else
              CCoMx = A(1)**ix
              CCoMxd = real(ix,kind=wp)*A(1)**(ix-1)
            end if

            if (iy == 0) then
              CCoMy = One
              CCoMyd = Zero
            else if (iy == 1) then
              CCoMy = A(2)
              CCoMyd = One
            else
              CCoMy = A(2)**iy
              CCoMyd = real(iy,kind=wp)*A(2)**(iy-1)
            end if

            if (iz == 0) then
              CCoMz = One
              CCoMzd = Zero
            else if (iz == 1) then
              CCoMz = A(3)
              CCoMzd = One
            else
              CCoMz = A(3)**iz
              CCoMzd = real(iz,kind=wp)*A(3)**(iz-1)
            end if
            tempd(1) = MM(ip,2)*ZA*CCoMxd*CCoMy*CCoMz
            tempd(2) = MM(ip,2)*ZA*CCoMx*CCoMyd*CCoMz
            tempd(3) = MM(ip,2)*ZA*CCoMx*CCoMy*CCoMzd
            if (iPrint >= 99) then
              write(u6,*) CCoMx,CCoMy,CCoMz
              write(u6,*) 'tempd=',tempd
            end if

            ! Distribute gradient

            nDisp = IndDsp(mdc+iCnt,iIrrep)
            do iCar=0,2
              iComp = 2**iCar
              if (TstFnc(dc(mdc+iCnt)%iCoSet,iIrrep,iComp,dc(mdc+iCnt)%nStab)) then
                nDisp = nDisp+1
                if (Direct(nDisp)) then
                  Temp(nDisp) = Temp(nDisp)-Tempd(iCar+1)
                end if
              end if
            end do

          end do
        end do

      end do
    end do
  end do
  if (iPrint >= 15) then
    Lab = ' The Nuclear Reaction Field (KirkWood) Contribution'
    call PrGrad(Lab,Temp,nGrad,ChDisp)
  end if

  call DaXpY_(nGrad,One,Temp,1,Grad,1)

else if (lRF .and. PCM) then

!***********************************************************************
!                                                                      *
!          Compute contributions due to the reaction field             *
!                      PCM Model                                       *
!                                                                      *
!***********************************************************************

  Temp(:) = Zero

  ! Loop over tiles

  do iTs=1,nTs
    ZA = PCM_SQ(1,iTs)+PCM_SQ(2,iTS)
    NoLoop = ZA == Zero
    ZA = ZA/real(nIrrep,kind=wp)
    if (NoLoop) cycle
    A(1:3) = PCMTess(1:3,iTs)

    ! Tile only stabilized by the unit operator

    nStb = 1
    iStb(0) = 0

    ndc = 0
    do jCnttp=1,nCnttp
      if (jCnttp > 1) ndc = ndc+dbsc(jCnttp-1)%nCntr
      ZB = dbsc(jCnttp)%Charge
      if (ZB == Zero) cycle
      if (dbsc(jCnttp)%pChrg) cycle
      if (dbsc(jCnttp)%Frag) cycle
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        PreFct = ZAZB*real(nIrrep,kind=wp)/real(LmbdR,kind=wp)
        do iR=0,nDCRR-1
          call OA(iDCRR(iR),B,RB)
          nOp = NrOpr(iDCRR(iR))
          if (EQ(A,RB)) cycle
          r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)

          ! The factor u/g will ensure that the value of the
          ! gradient in symmetry adapted and no symmetry basis
          ! will have the same value.

          fab = One
          dfab = Zero
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

          nDisp = IndDsp(ndc+jCnt,iIrrep)
          igv = nIrrep/dc(ndc+jCnt)%nStab
          do iCar=0,2
            dr_dB = -(A(iCar+1)-RB(iCar+1))/r12
            iComp = 2**iCar
            if (TstFnc(dc(ndc+jCnt)%iCoSet,iIrrep,iComp,dc(ndc+jCnt)%nStab)) then
              nDisp = nDisp+1
              if (Direct(nDisp)) then
                ps = real(iPrmt(nOp,iChBas(2+iCar)),kind=wp)
                Temp(nDisp) = Temp(nDisp)+ps*One/real(igv,kind=wp)*PreFct*dr_dB*df_dr
              end if
            end if
          end do   ! End loop over cartesian components, iCar

        end do     ! End loop over DCR operators, iR

      end do       ! End over centers, jCnt
    end do         ! End over basis set types, jCnttp
  end do           ! End of tiles

  if (iPrint >= 15) then
    Lab = ' The Nuclear Reaction Field (PCM) Contribution'
    call PrGrad(Lab,Temp,nGrad,ChDisp)
  end if

  call DaXpY_(nGrad,One,Temp,1,Grad,1)

  ! Add contribution due to the tiles.

  Temp(:) = Zero
  call PCM_Cav_grd(Temp,nGrad)
  if (iPrint >= 15) then
    Lab = ' The Cavity PCM Contribution'
    call PrGrad(Lab,Temp,nGrad,ChDisp)
  end if
  call DaXpY_(nGrad,One,Temp,1,Grad,1)

  ! Add contribution due to the electric field on the tiles.

  if (Conductor) then
    Temp(:) = Zero
    call PCM_EF_grd(Temp,nGrad)
    if (iPrint >= 15) then
      Lab = ' The EF PCM Contribution'
      call PrGrad(Lab,Temp,nGrad,ChDisp)
    end if
    call DaXpY_(nGrad,-One,Temp,1,Grad,1)
  end if

end if

return

end subroutine DrvN1
