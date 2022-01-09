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

use external_centers
use Basis_Info
use Center_Info
use Phase_Info
use Temporary_Parameters, only: Expert
use Real_Info, only: PotNuc
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 A(3), B(3), RB(3)
integer iDCRR(0:7), jCoSet(8,8), iStb(0:7), jStb(0:7)
logical EQ, NoLoop

NoLoop = .true.
iDum = 0
r12_Min = 0.0d0

! Nuclear repulsion, in case of some ECP we include the core electronic
! contribution (pseudo charges). The interaction of pseudo charges is
! excluded from the energy term.

PotNuc = Zero
mdc = 0
ZB = Zero
do iCnttp=1,nCnttp
  ZA = dbsc(iCnttp)%Charge
  if (dbsc(iCnttp)%Frag) ZA = dbsc(iCnttp)%FragCharge
  if (ZA == Zero) Go To 101
  do iCnt=1,dbsc(iCnttp)%nCntr
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    ndc = 0
    do jCnttp=1,iCnttp
      if (dbsc(iCnttp)%pChrg .and. dbsc(jCnttp)%pChrg) Go To 201
      if (dbsc(iCnttp)%Frag .and. dbsc(jCnttp)%Frag) Go To 201
      ZB = dbsc(jCnttp)%Charge
      if (dbsc(jCnttp)%Frag) ZB = dbsc(jCnttp)%FragCharge
      if (ZB == Zero) Go To 201
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
          RB(1) = dble(iPhase(1,iDCRR(iR)))*B(1)
          RB(2) = dble(iPhase(2,iDCRR(iR)))*B(2)
          RB(3) = dble(iPhase(3,iDCRR(iR)))*B(3)
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
                Gamma = dbsc(iCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gamma*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(iCnttp)%nM2
                Gamma = dbsc(iCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gamma*r12**2)
              end do
            end if
            if (dbsc(jCnttp)%ECP) then
              ! Add contribution from M1 operator
              do iM1xp=1,dbsc(jCnttp)%nM1
                Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gamma*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(jCnttp)%nM2
                Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gamma*r12**2)
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
        PotNuc = PotNuc+(Fact*ZAZB*temp*dble(nIrrep))/dble(LmbdR)

        jxyz = jxyz+3
      end do
201   continue
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
  end do
101 continue
  mdc = mdc+dbsc(iCnttp)%nCntr
end do

if (Show) then
  write(6,*)
  write(6,'(11X,A,F16.8,A)') ' Nuclear Potential Energy        ',PotNuc,' au'
  write(6,*)
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
    if (NoLoop) Go To 102
    A(1:3) = XF(1:3,iFd)
    iChxyz = iChAtm(A)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    ndc = 0
    do jCnttp=1,nCnttp
      ZB = dbsc(jCnttp)%Charge
      if (dbsc(jCnttp)%pChrg) Go To 202
      if (ZB == Zero) Go To 202
      if (dbsc(jCnttp)%Frag) Go To 202
      ZAZB = ZA*ZB
      do jCnt=1,dbsc(jCnttp)%nCntr
        B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

        ! Find the DCR for the two centers

        call DCR(LmbdR,iStb,nStb,dc(ndc+jCnt)%iStab,dc(ndc+jCnt)%nStab,iDCRR,nDCRR)

        temp0 = Zero
        temp1 = Zero
        temp2 = Zero
        do iR=0,nDCRR-1
          RB(1) = dble(iPhase(1,iDCRR(iR)))*B(1)
          RB(2) = dble(iPhase(2,iDCRR(iR)))*B(2)
          RB(3) = dble(iPhase(3,iDCRR(iR)))*B(3)
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
                Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                fab = fab+CffM1*exp(-Gamma*r12**2)
              end do
              ! Add contribution from M2 operator
              do iM2xp=1,dbsc(jCnttp)%nM2
                Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                fab = fab+CffM2*r12*exp(-Gamma*r12**2)
              end do
            end if
            temp0 = temp0+fab/r12
            if (nOrd_XF >= 1) temp1 = temp1-fab*(DAx*ABx+DAy*ABy+DAz*ABz)/r12**3
            if (nOrd_XF >= 2) then

              temp2 = temp2+fab*0.5d0*(3.0d0*(Qxx*ABx**2+2.0d0*Qxy*ABx*ABy+2.0d0*Qxz*ABx*ABz+Qyy*ABy**2+ &
                                              2.0d0*Qyz*ABy*ABz+Qzz*ABz**2)/r12**5-One/r12**3*(Qxx+Qyy+Qzz))
            end if

          end if
        end do
        PNX = PNX+((ZAZB*temp0+ZB*(temp1+temp2))*dble(nIrrep))/dble(LmbdR)

      end do
202   continue
      ndc = ndc+dbsc(jCnttp)%nCntr
    end do
102 continue
  end do

  if (Show) then
    write(6,'(19X,A,F16.8,A)') ' Nuclear-External Field Potential Energy ',PNX,' au'
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

  if (nXF > 9999) go to 1254

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
    elseif (nOrd_XF == 1) then
      ZA = XF(4,iFd)
      DAx = XF(5,iFd)
      DAy = XF(6,iFd)
      DAz = XF(7,iFd)
      NoLoop = (ZA == Zero) .and. (DAx == Zero) .and. (DAy == Zero) .and. (DAz == Zero)
    elseif (nOrd_XF == 2) then
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

    if (NoLoop) Go To 103
    A(1:3) = XF(1:3,iFd)
    iChxyz = iChAtm(A)
    call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)

    do jFd=1,iFd
      if (nOrd_XF == 0) then
        ZB = XF(4,jFd)
        NoLoop = ZB == Zero
      elseif (nOrd_XF == 1) then
        ZB = XF(4,jFd)
        DBx = XF(5,jFd)
        DBy = XF(6,jFd)
        DBz = XF(7,jFd)
        NoLoop = (ZB == Zero) .and. (DBx == Zero) .and. (DBy == Zero) .and. (DBz == Zero)
      elseif (nOrd_XF == 2) then
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

      if (NoLoop) Go To 203
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
        RB(1) = dble(iPhase(1,iDCRR(iR)))*B(1)
        RB(2) = dble(iPhase(2,iDCRR(iR)))*B(2)
        RB(3) = dble(iPhase(3,iDCRR(iR)))*B(3)
        DRBx = dble(iPhase(1,iDCRR(iR)))*DBx
        DRBy = dble(iPhase(2,iDCRR(iR)))*DBy
        DRBz = dble(iPhase(3,iDCRR(iR)))*DBz
        QRBxx = QBxx
        QRByy = QByy
        QRBzz = QBzz
        QRBxy = dble(iPhase(1,iDCRR(iR))*iPhase(2,iDCRR(iR)))*QBxy
        QRBxz = dble(iPhase(1,iDCRR(iR))*iPhase(3,iDCRR(iR)))*QBxz
        QRByz = dble(iPhase(2,iDCRR(iR))*iPhase(3,iDCRR(iR)))*QByz

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

            QAsum = (QAxx*x*x+QAyy*y*y+QAzz*z*z+2.0d0*(QAxy*x*y+QAxz*x*z+QAyz*y*z))
            QBsum = (QRBxx*x*x+QRByy*y*y+QRBzz*z*z+2.0d0*(QRBxy*x*y+QRBxz*x*z+QRByz*y*z))

            ! Q-Z
            eZQ = 0.5d0*3.0d0/r12**5*ZA*QBsum-0.5d0/r12**3*ZA*(QRBxx+QRByy+QRBzz)
            eQZ = 0.5d0*3.0d0/r12**5*ZB*QAsum-0.5d0/r12**3*ZB*(QAxx+QAyy+QAzz)

            ! Q-D
            eDQ = 0.5d0*(-15.0d0/r12**7*(DAx*x+DAy*y+DAz*z)*QBsum+3.0d0/r12**5* &
                  (3.0d0*DAx*QRBxx*x+DAy*QRBxx*y+DAz*QRBxx*z+2.0d0*DAx*QRBxy*y+2.0d0*DAy*QRBxy*x+2.0d0*DAx*QRBxz*z+ &
                   2.0d0*DAz*QRBxz*x+DAx*QRByy*x+3.0d0*DAy*QRByy*y+DAz*QRByy*z+2.0d0*DAy*QRByz*z+2.0d0*DAz*QRByz*y+DAx*QRBzz*x+ &
                   DAy*QRBzz*y+3.0d0*DAz*QRBzz*z))

            eQD = -0.5d0*(-15.0d0/r12**7*(DRBx*x+DRBy*y+DRBz*z)*QAsum+3.0d0/r12**5* &
                  (3.0d0*DRBx*QAxx*x+DRBy*QAxx*y+DRBz*QAxx*z+2.0d0*DRBx*QAxy*y+2.0d0*DRBy*QAxy*x+2.0d0*DRBx*QAxz*z+ &
                   2.0d0*DRBz*QAxz*x+DRBx*QAyy*x+3.0d0*DRBy*QAyy*y+DRBz*QAyy*z+2.0d0*DRBy*QAyz*z+2.0d0*DRBz*QAyz*y+DRBx*QAzz*x+ &
                   DRBy*QAzz*y+3.0d0*DRBz*QAzz*z))

            ! Q-Q
            eQQ = 0.25d0*(105.0d0/r12**9*QAsum*QBsum-15.0d0/r12**7* &
                  (6.0d0*QAxx*QRBxx*x*x+6.0d0*QAxx*QRBxy*x*y+6.0d0*QAxx*QRBxz*x*z+QAxx*QRByy*x*x+QAxx*QRByy*y*y+ &
                   2.0d0*QAxx*QRByz*y*z+QAxx*QRBzz*x*x+QAxx*QRBzz*z*z+6.0d0*QAxy*QRBxx*x*y+4.0d0*QAxy*QRBxy*x*x+ &
                   4.0d0*QAxy*QRBxy*y*y+4.0d0*QAxy*QRBxz*y*z+6.0d0*QAxy*QRByy*x*y+4.0d0*QAxy*QRByz*x*z+2.0d0*QAxy*QRBzz*x*y+ &
                   6.0d0*QAxz*QRBxx*x*z+4.0d0*QAxz*QRBxy*y*z+4.0d0*QAxz*QRBxz*x*x+4.0d0*QAxz*QRBxz*z*z+2.0d0*QAxz*QRByy*x*z+ &
                   4.0d0*QAxz*QRByz*x*y+6.0d0*QAxz*QRBzz*x*z+QAyy*QRBxx*x*x+QAyy*QRBxx*y*y+6.0d0*QAyy*QRBxy*x*y+ &
                   2.0d0*QAyy*QRBxz*x*z+6.0d0*QAyy*QRByy*y*y+6.0d0*QAyy*QRByz*y*z+QAyy*QRBzz*y*y+QAyy*QRBzz*z*z+ &
                   2.0d0*QAyz*QRBxx*y*z+4.0d0*QAyz*QRBxy*x*z+4.0d0*QAyz*QRBxz*x*y+6.0d0*QAyz*QRByy*y*z+4.0d0*QAyz*QRByz*y*y+ &
                   4.0d0*QAyz*QRByz*z*z+6.0d0*QAyz*QRBzz*y*z+QAzz*QRBxx*x*x+QAzz*QRBxx*z*z+2.0d0*QAzz*QRBxy*x*y+ &
                   6.0d0*QAzz*QRBxz*x*z+QAzz*QRByy*y*y+QAzz*QRByy*z*z+6.0d0*QAzz*QRByz*y*z+6.0d0*QAzz*QRBzz*z*z)+ &
                  3.0d0/r12**5* &
                  (3.0d0*QAxx*QRBxx+QAxx*QRByy+QAxx*QRBzz+4.0d0*QAxy*QRBxy+4.0d0*QAxz*QRBxz+QAyy*QRBxx+ &
                   3.0d0*QAyy*QRByy+QAyy*QRBzz+4.0d0*QAyz*QRByz+QAzz*QRBxx+QAzz*QRByy+3.0d0*QAzz*QRBzz))

            eTOT = eZZ+eZD+eDZ+eDD+eZQ+eQZ+eDQ+eQD+eQQ
            !write(6,*)'eZZ',eZZ
            !write(6,*)'eDZ',eDZ
            !write(6,*)'eZD',eZD
            !write(6,*)'eDD',eDD
            !write(6,*)'eZQ',eZQ
            !write(6,*)'eQZ',eQZ
            !write(6,*)'eDQ',eDQ
            !write(6,*)'eQD',eQD
            !write(6,*)'eQQ',eQQ
            !write(6,*)'eTOT',eTOT

          end if ! if quadrupoles are present

          temp = temp+eTOT

        end if
      end do
      !PXX = PXX+(Fact*(ZAZB*temp0+ZB*temp1+ZA*temp2+temp3)*dble(nIrrep))/dble(LmbdR)
      PXX = PXX+(Fact*temp*dble(nIrrep))/dble(LmbdR)

203   continue
    end do
103 continue
  end do

  if (Show) then
    write(6,'(19X,A,F16.8,A)') ' External Field Potential Energy         ',PXX,' au '
  end if
  call Put_dScalar('PC Self Energy',PXX)

  !PotNuc = PotNuc+PXX

  if (Show) then
    write(6,*)
    write(6,*)
    write(6,'(11X,A,F16.8,A)') ' Total Nuclear Potential Energy        ',PotNuc,' au'
    write(6,*)
  end if

1254 continue
end if

call Put_dScalar('PotNuc',PotNuc)
if (isstructure() == 1) then
  call Add_Info('PotNuc',[PotNuc],1,6)
else
  call Add_Info('PotNuc',[PotNuc],1,12)
end if

return

end subroutine DrvN0
