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

subroutine DrvN0()
!***********************************************************************
!                                                                      *
! Object: to compute the nuclear contributions to the nuclear potential*
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!             Modified for various other contributions May 95', RL     *
!***********************************************************************

use External_Centers, only: XF, nOrd_XF, nXF
use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Phase_Info, only: iPhase
use Gateway_global, only: Expert
use Gateway_Info, only: PotNuc
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One, Two, Three, Four, Six, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iChxyz, iCnt, iCnttp, iDCRR(0:7), iDum, iFd, iM1xp, iM2xp, iR, iStb(0:7), jCnt, jCntMx, jCnttp, jCoSet(8,8), &
                     jFd, jStb(0:7), LmbdR, mdc, mStb, ndc, nDCRR, nStb
real(kind=wp) :: A(3), ABx, ABy, ABz, B(3), CffM1, CffM2, DAx, DAy, DAz, DBx, DBy, DBz, DRBx, DRBy, DRBz, eDD, eDQ, eDZ, eQD, eQQ, &
                 eQZ, eTot, eZD, eZQ, eZZ, fab, Fact, Gam, PNX, PXX, QAsum, QAxx, QAxy, QAxz, QAyy, QAyz, QAzz, QBsum, QBxx, QBxy, &
                 QBxz, QByy, QByz, QBzz, QRBxx, QRBxy, QRBxz, QRByy, QRByz, QRBzz, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, r12, r12_Min, &
                 RB(3), temp, temp0, temp1, temp2, x, y, z, ZA, ZAZB, ZB
logical(kind=iwp) :: EQ, NoLoop
integer(kind=iwp), external :: iChAtm, isstructure

NoLoop = .true.
iDum = 0
r12_Min = Zero

! Nuclear repulsion, in case of some ECP we include the core electronic
! contribution (pseudo charges). The interaction of pseudo charges is
! excluded from the energy term.

PotNuc = Zero
mdc = 0
ZB = Zero
do iCnttp=1,nCnttp
  ZA = dbsc(iCnttp)%Charge
  if (dbsc(iCnttp)%Frag) ZA = dbsc(iCnttp)%FragCharge
  mdc = mdc+dbsc(iCnttp)%nCntr
  if (ZA == Zero) cycle
  mdc = mdc-dbsc(iCnttp)%nCntr
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=1,iCnttp
      ndc = ndc+dbsc(jCnttp)%nCntr
      if (dbsc(iCnttp)%pChrg .and. dbsc(jCnttp)%pChrg) cycle
      if (dbsc(iCnttp)%Frag .and. dbsc(jCnttp)%Frag) cycle
      ZB = dbsc(jCnttp)%Charge
      if (dbsc(jCnttp)%Frag) ZB = dbsc(jCnttp)%FragCharge
      if (ZB == Zero) cycle
      ndc = ndc-dbsc(jCnttp)%nCntr
      ZAZB = ZA*ZB
      jCntMx = dbsc(jCnttp)%nCntr
      if (iCnttp == jCnttp) jCntMx = iCnt
      do jCnt=1,jCntMx
        ! Introduce factor to ensure that contributions from
        ! A>B are the only to be accumulated.
        Fact = One
        if ((iCnttp == jCnttp) .and. (iCnt == jCnt)) Fact = Half
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        temp = Zero
        do iR=0,nDCRR-1
          RB(1) = real(iPhase(1,iDCRR(iR)),kind=wp)*B(1)
          RB(2) = real(iPhase(2,iDCRR(iR)),kind=wp)*B(2)
          RB(3) = real(iPhase(3,iDCRR(iR)),kind=wp)*B(3)
          ! The index A=RB is illegal.
          if (.not. EQ(A,RB)) then
            r12 = sqrt((A(1)-RB(1))**2+(A(2)-RB(2))**2+(A(3)-RB(3))**2)
            if ((r12 < r12_min) .and. (.not. Expert)) then
              call WarningMessage(2,' The distance between two centers are found to be unphysically too short.; '// &
                                  'If you want to persist operate in expert mode.')
              call Abend()
            end if
            fab = One
            if (dbsc(iCnttp)%ECP) then
              ! Add contribution from M1 operator
              do iM1xp=1,dbsc(iCnttp)%nM1
                Gam = dbsc(iCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gam*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(iCnttp)%nM2
                Gam = dbsc(iCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gam*r12**2)
              end do
            end if
            if (dbsc(jCnttp)%ECP) then
              ! Add contribution from M1 operator
              do iM1xp=1,dbsc(jCnttp)%nM1
                Gam = dbsc(jCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gam*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(jCnttp)%nM2
                Gam = dbsc(jCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gam*r12**2)
              end do
            end if

            temp = temp+fab/r12
          else
            ! turn off checking if in expert mode, for wavefunction overlaps
            if (.not. Expert) then
              if ((iCnttp /= jCnttp) .or. (iCnt /= jCnt)) then
                call WarningMessage(2,'You have two charges on top of each other!; Correct the input!')
                call Quit_OnUserError()
              end if
            end if
          end if
        end do
        PotNuc = PotNuc+(Fact*ZAZB*temp*real(nIrrep,kind=wp))/real(LmbdR,kind=wp)
      end do
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
  end do
  mdc = mdc+dbsc(iCnttp)%nCntr
end do

if (Show) then
  write(u6,*)
  write(u6,'(11X,A,F16.8,A)') ' Nuclear Potential Energy        ',PotNuc,' au'
  write(u6,*)
end if

if (allocated(XF) .and. (nOrd_XF >= 0)) then

  ! Add contribution for interaction external field and nuclear
  ! charges. Here we will have charge-charge, and charge-dipole
  ! inteaction.

  ZA = Zero
  DAx = Zero
  DAy = Zero
  DAz = Zero
  Qxx = Zero
  Qxy = Zero
  Qxz = Zero
  Qyy = Zero
  Qyz = Zero
  Qzz = Zero

  PNX = Zero
  iDum = 0
  do iFd=1,nXF
    if (nOrd_XF == 0) then
      ZA = XF(4,iFd)
      NoLoop = ZA == Zero
    else if (nOrd_XF == 1) then
      ZA = XF(4,iFd)
      DAx = XF(5,iFd)
      DAy = XF(6,iFd)
      DAz = XF(7,iFd)
      NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero)
    else if (nOrd_XF == 2) then
      ZA = XF(4,iFd)
      DAx = XF(5,iFd)
      DAy = XF(6,iFd)
      DAz = XF(7,iFd)
      Qxx = XF(8,iFd)
      Qxy = XF(9,iFd)
      Qxz = XF(10,iFd)
      Qyy = XF(11,iFd)
      Qyz = XF(12,iFd)
      Qzz = XF(13,iFd)
      NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero) .and. (Qxx == Zero) .and. &
               (Qxy == Zero) .and. (Qxz == Zero) .and. (Qyy == Zero) .and. (Qyz == Zero) .and. (Qzz == Zero)
    else
      call WarningMessage(2,'Option not implemented yet!')
      call Quit_OnUserError()
    end if
    if (NoLoop) cycle
    A(1:3) = XF(1:3,iFd)
    iChxyz = iChAtm(A)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    ndc = 0
    do jCnttp=1,nCnttp
      ZB = dbsc(jCnttp)%Charge
      ndc = ndc+dbsc(jCnttp)%nCntr
      if (dbsc(jCnttp)%pChrg) cycle
      if (ZB == Zero) cycle
      if (dbsc(jCnttp)%Frag) cycle
      ndc = ndc-dbsc(jCnttp)%nCntr
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        temp0 = Zero
        temp1 = Zero
        temp2 = Zero
        do iR=0,nDCRR-1
          RB(1) = real(iPhase(1,iDCRR(iR)),kind=wp)*B(1)
          RB(2) = real(iPhase(2,iDCRR(iR)),kind=wp)*B(2)
          RB(3) = real(iPhase(3,iDCRR(iR)),kind=wp)*B(3)
          ! The index A=RB is illegal.
          if (.not. EQ(A,RB)) then
            ABx = A(1)-RB(1)
            ABy = A(2)-RB(2)
            ABz = A(3)-RB(3)
            r12 = sqrt(ABx**2+ABy**2+ABz**2)

            fab = One
            if (dbsc(jCnttp)%ECP) then
              ! Add contribution from M1 operator
              do iM1xp=1,dbsc(jCnttp)%nM1
                Gam = dbsc(jCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gam*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(jCnttp)%nM2
                Gam = dbsc(jCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gam*r12**2)
              end do
            end if
            temp0 = temp0+fab/r12
            if (nOrd_XF >= 1) temp1 = temp1-fab*(DAx*ABx+DAy*ABy+DAz*ABz)/r12**3
            if (nOrd_XF >= 2) then

              temp2 = temp2+fab*Half*(Three*(Qxx*ABx**2+Two*Qxy*ABx*ABy+Two*Qxz*ABx*ABz+Qyy*ABy**2+ &
                                             Two*Qyz*ABy*ABz+Qzz*ABz**2)/r12**5-One/r12**3*(Qxx+Qyy+Qzz))
            end if

          end if
        end do
        PNX = PNX+((ZAZB*temp0+ZB*(temp1+temp2))*real(nIrrep,kind=wp))/real(LmbdR,kind=wp)

      end do
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
  end do

  if (Show) then
    write(u6,'(19X,A,F16.8,A)') ' Nuclear-External Field Potential Energy ',PNX,' au'
  end if

  PotNuc = PotNuc+PNX

  ! Add contribution for self interaction of the external field.
  ! Here we will have charge-charge, charge-dipole, and dipole-
  ! dipole inteaction. This term will only be computed but not
  ! added to the acutual nuclear potential energy.

  ! LU: in the case when the number of point charges is larger than
  !     10^4, the algorith below is extremely slow, as it implies a
  !     double loop over all field points (nXF). The resulting value
  !     is stored in the RunFile (label='PC Self Energy'), but never
  !     used again in the entire MOLCAS code.
  !     As a first attempt to optimize this part of the code, we skip
  !     the computation of the self interaction of the external field.

  if (nXF <= 9999) then

    PXX = Zero

    DAx = Zero
    DAy = Zero
    DAz = Zero
    QAxx = Zero
    QAxy = Zero
    QAxz = Zero
    QAyy = Zero
    QAyz = Zero
    QAzz = Zero
    DBx = Zero
    DBy = Zero
    DBz = Zero
    QBxx = Zero
    QBxy = Zero
    QBxz = Zero
    QByy = Zero
    QByz = Zero
    QBzz = Zero

    iDum = 0
    do iFd=1,nXF
      if (nOrd_XF == 0) then
        ZA = XF(4,iFd)
        NoLoop = ZA == Zero
      else if (nOrd_XF == 1) then
        ZA = XF(4,iFd)
        DAx = XF(5,iFd)
        DAy = XF(6,iFd)
        DAz = XF(7,iFd)
        NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero)
      else if (nOrd_XF == 2) then
        ZA = XF(4,iFd)
        DAx = XF(5,iFd)
        DAy = XF(6,iFd)
        DAz = XF(7,iFd)
        QAxx = XF(8,iFd)
        QAxy = XF(9,iFd)
        QAxz = XF(10,iFd)
        QAyy = XF(11,iFd)
        QAyz = XF(12,iFd)
        QAzz = XF(13,iFd)
        NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero) .and. (QAxx == Zero) .and. &
                 (QAxy == Zero) .and. (QAxz == Zero) .and. (QAyy == Zero) .and. (QAyz == Zero) .and. (QAzz == Zero)
      else
        call WarningMessage(2,'Option not implemented yet!')
        call Quit_OnUserError()
      end if

      if (NoLoop) cycle
      A(1:3) = XF(1:3,iFd)
      iChxyz = iChAtm(A)
      call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

      do jFd=1,iFd
        if (nOrd_XF == 0) then
          ZB = XF(4,jFd)
          NoLoop = ZB == Zero
        else if (nOrd_XF == 1) then
          ZB = XF(4,jFd)
          DBx = XF(5,jFd)
          DBy = XF(6,jFd)
          DBz = XF(7,jFd)
          NoLoop = (ZB == Zero) .and. (DBx == Zero) .and. (DBy == Zero) .and. (DBz == Zero)
        else if (nOrd_XF == 2) then
          ZB = XF(4,jFd)
          DBx = XF(5,jFd)
          DBy = XF(6,jFd)
          DBz = XF(7,jFd)
          QBxx = XF(8,jFd)
          QBxy = XF(9,jFd)
          QBxz = XF(10,jFd)
          QByy = XF(11,jFd)
          QByz = XF(12,jFd)
          QBzz = XF(13,jFd)
          NoLoop = (ZB == Zero) .and. (DBx == Zero) .and. (DBy == Zero) .and. (DBz == Zero) .and. (QBxx == Zero) .and. &
                   (QBxy == Zero) .and. (QBxz == Zero) .and. (QByy == Zero) .and. (QByz == Zero) .and. (QBzz == Zero)
        else
          call WarningMessage(2,'Option not implemented yet!')
          call Quit_OnUserError()
        end if

        if (NoLoop) cycle
        ZAZB = ZA*ZB
        B(1:3) = XF(1:3,jFd)
        iChxyz = iChAtm(B)
        call Stblz(iChxyz,mStb,jStb,iDum,jCoSet)
        ! Introduce factor to ensure that contributions from
        ! A>B are the only to be accumulated.
        Fact = One
        if (iFd == jFd) Fact = Half

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,jStb,mStb,iDCRR,nDCRR)

        temp = Zero
        do iR=0,nDCRR-1
          RB(1) = real(iPhase(1,iDCRR(iR)),kind=wp)*B(1)
          RB(2) = real(iPhase(2,iDCRR(iR)),kind=wp)*B(2)
          RB(3) = real(iPhase(3,iDCRR(iR)),kind=wp)*B(3)
          DRBx = real(iPhase(1,iDCRR(iR)),kind=wp)*DBx
          DRBy = real(iPhase(2,iDCRR(iR)),kind=wp)*DBy
          DRBz = real(iPhase(3,iDCRR(iR)),kind=wp)*DBz
          QRBxx = QBxx
          QRByy = QByy
          QRBzz = QBzz
          QRBxy = real(iPhase(1,iDCRR(iR))*iPhase(2,iDCRR(iR)),kind=wp)*QBxy
          QRBxz = real(iPhase(1,iDCRR(iR))*iPhase(3,iDCRR(iR)),kind=wp)*QBxz
          QRByz = real(iPhase(2,iDCRR(iR))*iPhase(3,iDCRR(iR)),kind=wp)*QByz

          ! The index A=RB is illegal.
          if (.not. EQ(A,RB)) then
            x = A(1)-RB(1)
            y = A(2)-RB(2)
            z = A(3)-RB(3)
            r12 = sqrt(x**2+y**2+z**2)

            eZZ = ZAZB*One/r12
            eDZ = ZB*(-(DAx*x+DAy*y+DAz*z))/r12**3
            eZD = ZA*(DRBx*x+DRBy*y+DRBz*z)/r12**3
            eDD = (DAx*DRBx+DAy*DRBy+DAz*DRBz)/r12**3-Three*(DAx*x+DAy*y+DAz*z)*(DRBx*x+DRBy*y+DRBz*z)/r12**5

            if (nOrd_XF < 2) then
              eTOT = eZZ+eZD+eDZ+eDD
            else

              QAsum = (QAxx*x*x+QAyy*y*y+QAzz*z*z+Two*(QAxy*x*y+QAxz*x*z+QAyz*y*z))
              QBsum = (QRBxx*x*x+QRByy*y*y+QRBzz*z*z+Two*(QRBxy*x*y+QRBxz*x*z+QRByz*y*z))

              ! Q-Z
              eZQ = Half*Three/r12**5*ZA*QBsum-Half/r12**3*ZA*(QRBxx+QRByy+QRBzz)
              eQZ = Half*Three/r12**5*ZB*QAsum-Half/r12**3*ZB*(QAxx+QAyy+QAzz)

              ! Q-D
              eDQ = Half*(-15.0_wp/r12**7*(DAx*x+DAy*y+DAz*z)*QBsum+Three/r12**5* &
                    (Three*DAx*QRBxx*x+DAy*QRBxx*y+DAz*QRBxx*z+Two*DAx*QRBxy*y+Two*DAy*QRBxy*x+Two*DAx*QRBxz*z+Two*DAz*QRBxz*x+ &
                     DAx*QRByy*x+Three*DAy*QRByy*y+DAz*QRByy*z+Two*DAy*QRByz*z+Two*DAz*QRByz*y+DAx*QRBzz*x+DAy*QRBzz*y+ &
                     Three*DAz*QRBzz*z))

              eQD = -Half*(-15.0_wp/r12**7*(DRBx*x+DRBy*y+DRBz*z)*QAsum+Three/r12**5* &
                    (Three*DRBx*QAxx*x+DRBy*QAxx*y+DRBz*QAxx*z+Two*DRBx*QAxy*y+Two*DRBy*QAxy*x+Two*DRBx*QAxz*z+Two*DRBz*QAxz*x+ &
                     DRBx*QAyy*x+Three*DRBy*QAyy*y+DRBz*QAyy*z+Two*DRBy*QAyz*z+Two*DRBz*QAyz*y+DRBx*QAzz*x+DRBy*QAzz*y+ &
                     Three*DRBz*QAzz*z))

              ! Q-Q
              eQQ = Quart*(105.0_wp/r12**9*QAsum*QBsum-15.0_wp/r12**7* &
                    (Six*QAxx*QRBxx*x*x+Six*QAxx*QRBxy*x*y+Six*QAxx*QRBxz*x*z+QAxx*QRByy*x*x+QAxx*QRByy*y*y+Two*QAxx*QRByz*y*z+ &
                     QAxx*QRBzz*x*x+QAxx*QRBzz*z*z+Six*QAxy*QRBxx*x*y+Four*QAxy*QRBxy*x*x+Four*QAxy*QRBxy*y*y+Four*QAxy*QRBxz*y*z+ &
                     Six*QAxy*QRByy*x*y+Four*QAxy*QRByz*x*z+Two*QAxy*QRBzz*x*y+Six*QAxz*QRBxx*x*z+Four*QAxz*QRBxy*y*z+ &
                     Four*QAxz*QRBxz*x*x+Four*QAxz*QRBxz*z*z+Two*QAxz*QRByy*x*z+Four*QAxz*QRByz*x*y+Six*QAxz*QRBzz*x*z+ &
                     QAyy*QRBxx*x*x+QAyy*QRBxx*y*y+Six*QAyy*QRBxy*x*y+Two*QAyy*QRBxz*x*z+Six*QAyy*QRByy*y*y+Six*QAyy*QRByz*y*z+ &
                     QAyy*QRBzz*y*y+QAyy*QRBzz*z*z+Two*QAyz*QRBxx*y*z+Four*QAyz*QRBxy*x*z+Four*QAyz*QRBxz*x*y+Six*QAyz*QRByy*y*z+ &
                     Four*QAyz*QRByz*y*y+Four*QAyz*QRByz*z*z+Six*QAyz*QRBzz*y*z+QAzz*QRBxx*x*x+QAzz*QRBxx*z*z+Two*QAzz*QRBxy*x*y+ &
                     Six*QAzz*QRBxz*x*z+QAzz*QRByy*y*y+QAzz*QRByy*z*z+Six*QAzz*QRByz*y*z+Six*QAzz*QRBzz*z*z)+ &
                    Three/r12**5* &
                    (Three*QAxx*QRBxx+QAxx*QRByy+QAxx*QRBzz+Four*QAxy*QRBxy+Four*QAxz*QRBxz+QAyy*QRBxx+ &
                     Three*QAyy*QRByy+QAyy*QRBzz+Four*QAyz*QRByz+QAzz*QRBxx+QAzz*QRByy+Three*QAzz*QRBzz))

              eTOT = eZZ+eZD+eDZ+eDD+eZQ+eQZ+eDQ+eQD+eQQ
              !write(u6,*)'eZZ',eZZ
              !write(u6,*)'eDZ',eDZ
              !write(u6,*)'eZD',eZD
              !write(u6,*)'eDD',eDD
              !write(u6,*)'eZQ',eZQ
              !write(u6,*)'eQZ',eQZ
              !write(u6,*)'eDQ',eDQ
              !write(u6,*)'eQD',eQD
              !write(u6,*)'eQQ',eQQ
              !write(u6,*)'eTOT',eTOT

            end if ! if quadrupoles are present

            temp = temp+eTOT

          end if
        end do
        !PXX = PXX+(Fact*(ZAZB*temp0+ZB*temp1+ZA*temp2+temp3)*real(nIrrep,kind=wp))/real(LmbdR,kind=wp)
        PXX = PXX+(Fact*temp*real(nIrrep,kind=wp))/real(LmbdR,kind=wp)

      end do
    end do

    if (Show) then
      write(u6,'(19X,A,F16.8,A)') ' External Field Potential Energy         ',PXX,' au '
    end if
    call Put_dScalar('PC Self Energy',PXX)

    !PotNuc = PotNuc+PXX

    if (Show) then
      write(u6,*)
      write(u6,*)
      write(u6,'(11X,A,F16.8,A)') ' Total Nuclear Potential Energy        ',PotNuc,' au'
      write(u6,*)
    end if

  end if
end if

call Put_dScalar('PotNuc',PotNuc)
if (isstructure() == 1) then
  call Add_Info('PotNuc',[PotNuc],1,6)
else
  call Add_Info('PotNuc',[PotNuc],1,12)
end if

return

end subroutine DrvN0
