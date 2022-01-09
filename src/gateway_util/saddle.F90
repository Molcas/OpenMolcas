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
! Copyright (C) 2009, Roland Lindh                                     *
!               2010, Mickael G. Delcey                                *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Saddle()
!***********************************************************************
!                                                                      *
! Object: to set up for a TS optimization with the SADDLE approach.    *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             January 2009                                             *
!***********************************************************************

use Basis_Info
use Center_Info
use external_centers
use Sizes_of_Seward, only: S
use Real_Info, only: E1, E2, SadStep, Shake
use Logical_Info, only: Align_Only, Do_Align, lRP, lRP_Post
use Symmetry_Info, only: nIrrep, VarR, VarT

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
character*1 Mode
character*16 StdIn
logical Not_First_Iter, FindTS, Found, quadratic, Invar
#include "angstr.fh"
#include "warnings.h"
character*2, dimension(:), allocatable :: Elm
real*8, dimension(:), allocatable :: TanVec, TmpA, W
real*8, dimension(:,:), allocatable :: Vec, MEP
integer, dimension(:), allocatable :: iStab
integer ipX2, ipX3
integer iRef, iOpt
real*8, dimension(:,:), allocatable :: XYZ
real*8 RandVect(3)
#include "periodic_table.fh"

!***********************************************************************
!                                                                      *
!                            Prologue                                  *
!                                                                      *
!***********************************************************************

nSaddle_Max = 100
quadratic = .false.

! Avoid warnings

R11 = Zero
R22 = Zero
R1R2 = Zero
iX0 = -1
iX1 = -1

! If lRP true in Info and Saddle block active but set to zero then
! this is after the Saddle procedure is terminated and lRP should
! be ignored.

if (lRP) then
  call qpg_dArray('Saddle',Not_First_Iter,nData)
  if (Not_First_Iter .and. (nData == 0)) lRP = .false.
end if
lRP_Post = lRP
!                                                                      *
!***********************************************************************
!                                                                      *
! Get informations about Saddle in RunFile

if (lRP) then

  ! Find out whether the energy is assumed invariant to trans. & rot.

  Invar = .not. (VarR .or. VarT)

  nAt = nRP/3
  call qpg_dArray('Saddle',Not_First_Iter,nData)
  nSaddle = 2*nRP+5
  call mma_allocate(TmpA,nSaddle,label='TmpA')
  if (Not_First_Iter) then
    call Get_dArray('Saddle',TmpA,nSaddle)
    Update = TmpA(2*nRP+5)

    if (Update /= One) then
      HSR = TmpA(2*nRP+4)
      call mma_deallocate(TmpA)
      lRP_Post = .false.
      Go To 100
    end if
    Update = Zero

    call dcopy_(3*nAt,TmpA(1),1,RP_Centers(1,1,1),1)
    call dcopy_(3*nAt,TmpA(1+3*nAt),1,RP_Centers(1,1,2),1)
    E1 = TmpA(6*nAt+1)
    E2 = TmpA(6*nAt+2)
    HSR0 = TmpA(6*nAt+3)
  else
    dHSR = Zero ! Dummy initialize
    Update = Zero
    HSR0 = Zero ! Dummy initialize
    iOff_Iter = 0
    call Put_iScalar('iOff_Iter',iOff_Iter)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !       Get the symmetry stabilizers for each center

    call mma_allocate(iStab,nAt,label='iStab')
    iAt = 1
    nsc = 0
    do i=1,nCnttp
      do iCnt=1,dbsc(i)%nCntr
        nsc = nsc+1
        if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
          jTmp = 0
          do j=1,dc(nsc)%nStab-1
            jTmp = ior(jTmp,dc(nsc)%iStab(j))
          end do
          iStab(iAt) = jTmp
          iAt = iAt+1
        end if
      end do
    end do

    ! Shake initial structures

    if (Shake > Zero) then
      do iAt=1,nAt
        S%nDim = 0
        do j=0,2
          if (iand(iStab(iAt),2**j) == 0) S%nDim = S%nDim+1
        end do
        if (S%nDim > 0) then
          do iRP=1,2
            call Random_Vector(S%nDim,RandVect(1:S%nDim),.false.)
            jDim = 0
            do j=0,2
              if (iand(iStab(iAt),2**j) == 0) then
                jDim = jDim+1
                RP_Centers(jDim,iAt,iRP) = RP_Centers(jDim,iAt,iRP)+Shake*RandVect(jDim)
              end if
            end do
          end do
        end if
      end do
    end if

    ! Retrieve the weights even if the structures are not
    ! going to be aligned explicitly

    call mma_Allocate(XYZ,3*nAt*8,2,label='XYZ')
    iReac = 1
    iProd = 2
    call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iReac),mAt)
    call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iProd),mAt)
    call Qpg_dArray('Weights',Found,nData)
    if (Found .and. (nData >= mAt)) then
      call mma_allocate(W,nData,label='W')
      call Get_dArray('Weights',W,nData)
    else
      call SysAbendMsg('Saddle','No or wrong weights were found in the RUNFILE.','')
    end if

    ! Align the reactant and product structures the first time.
    ! Only if energy is invariant

    if (Do_Align .and. Invar) then
      ! Note: this might break symmetry
      call Superpose_w(XYZ(1,iReac),XYZ(1,iProd),W,mAt,RMS,RMSMax)
      call Fix_Symmetry(XYZ(1,iReac),nAt,iStab)
      call Add_Info('RMSD',[RMS],1,6)
      call Add_Info('RMSMax',[RMSMax],1,6)
      call dcopy_(3*nAt,XYZ(1,iReac),1,RP_Centers(1,1,1),1)
      call dcopy_(3*nAt,XYZ(1,iProd),1,RP_Centers(1,1,2),1)

      if (Align_Only) then

        ! Get the atomic symbols with symmetry unfolded

        call mma_allocate(Elm,nAt)
        iAt = 0
        iAtSym = nAt
        ndc = 0
        do iCnttp=1,nCnttp
          do iCnt=1,dbsc(iCnttp)%nCntr
            ndc = ndc+1
            if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
              iAt = iAt+1
              Elm(iAt) = PTab(dbsc(iCnttp)%AtmNr)
              do i=1,nIrrep/dc(ndc)%nStab-1
                iAtSym = iAtSym+1
                Elm(iAtSym) = Elm(iAt)
              end do
            end if
          end do
        end do

        write(6,*)
        write(6,*) 'Aligned Reactants and Products'
        write(6,*) '=============================='
        write(6,*)
        write(6,*)
        write(6,*) ' Reactants / Angstrom'
        write(6,*) '====================='
        write(6,*)
        do iAt=1,mAt
          write(6,'(A,1X,3F15.8)') Elm(iAt),(XYZ((iAt-1)*3+jAt,iReac)*Angstr,jAt=1,3)
        end do
        write(6,*)
        write(6,*)
        write(6,*) ' Products / Angstrom'
        write(6,*) '===================='
        write(6,*)
        do iAt=1,mAt
          write(6,'(A,1X,3F15.8)') Elm(iAt),(XYZ((iAt-1)*3+jAt,iProd)*Angstr,jAt=1,3)
        end do
        write(6,*)
        write(6,*)
        call WarningMessage(2,'Molecular alignment completed')
        write(6,*)
        iReturn = _RC_ALL_IS_WELL_
        call mma_deallocate(XYZ)
        call mma_deallocate(iStab)
        call mma_deallocate(W)
        call mma_deallocate(Elm)
        call ClsSew()
        call xQuit(iReturn)
      end if
    end if

    call mma_deallocate(XYZ)
    call mma_deallocate(iStab)
    call mma_deallocate(W)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Determine the distance between the two structures

  call mma_allocate(iStab,nAt,label='iStab')
  iAt = 1
  nsc = 0
  do i=1,nCnttp
    do iCnt=1,dbsc(i)%nCntr
      nsc = nsc+1
      if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
        jTmp = 0
        do j=1,dc(nsc)%nStab-1
          jTmp = ior(jTmp,dc(nsc)%iStab(j))
        end do
        iStab(iAt) = jTmp
        iAt = iAt+1
      end if
    end do
  end do

  call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
  iRA1 = 1
  iRA2 = 2
  call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt)
  call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt)
  call Qpg_dArray('Weights',Found,nData)
  if (Found .and. (nData >= mAt)) then
    call mma_allocate(W,nData,label='W')
    call Get_dArray('Weights',W,nData)
  else
    call SysAbendMsg('Saddle','No or wrong weights were found in the RUNFILE.','')
  end if
  if (.not. Invar) then

    ! If the energy is not trans/rot invariant, compute the weighted
    ! RMS with the current structures

    HSR = Zero
    wTot = Zero
    iOff = 1
    do i=1,mAt
      do ixyz=1,3
        diff = XYZ(iOff,iRA2)-XYZ(iOff,iRA1)
        HSR = HSR+W(i)*diff**2
        iOff = iOff+1
      end do
      wTot = wTot+W(i)
    end do
    HSR = sqrt(HSR/wTot)
  else
    call Get_RMSD_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,HSR)
  end if
  call mma_deallocate(XYZ)

  if (E1 <= E2) then
    Mode = 'R'
  else
    Mode = 'P'
  end if

  !*********************************************************************
  !                                                                    *
  !                 Determine the desired HSR                          *
  !                                                                    *
  !*********************************************************************
  if (Not_First_Iter) then

    ! FindTS

    FindTS = HSR <= (1.5d0*SadStep)
    if (FindTS) then
      write(6,*) '**************************'
      write(6,*) '* Enable TS optimization *'
      write(6,*) '**************************'
      Update = 2.0d0
      Delta = 0.5d0-6.25d0*(E2-E1)
      Delta = min(Delta,0.75d0)
      Delta = max(Delta,0.25d0)
      if (Mode == 'R') then
        HSR = (One-Delta)*HSR
      else
        HSR = Delta*HSR
      end if

      ! Calculate tangent vector

      call mma_allocate(TanVec,3*nAt,label='TanVec')
      if (Mode == 'R') then
        call Calc_LSTvec(3*nAt,RP_Centers(1,1,2),RP_Centers(1,1,1),TanVec,Invar)
      else
        call Calc_LSTvec(3*nAt,RP_Centers(1,1,1),RP_Centers(1,1,2),TanVec,Invar)
      end if
      call Put_dArray('TanVec',TanVec,3*nAt)
      call mma_deallocate(TanVec)
    else

      ! Try to be smart:
      ! be slow the first iterations and the last ones
      ! and also if large reorganisation the previous iteration

      call Qpg_iScalar('nMEP',Found)
      if (Found) then
        quadratic = .true.
        call Get_iScalar('nMEP',iSaddle)

        ! Compute some distances, used for quadratic interpolation

        call mma_Allocate(MEP,nRP,nSaddle_Max,label='MEP')
        call Get_dArray('MEP-Coor    ',MEP,nRP*nSaddle_Max)
        call mma_allocate(Vec,nRP,2,label='Vec')
        iX0 = iSaddle-2
        iX1 = iSaddle-1
        if (Mode == 'R') then
          ipX2 = 1
          ipX3 = 2
        else
          ipX2 = 2
          ipX3 = 1
        end if
        if (iSaddle < 3) iX0 = iX1

        ! Align everything with the current structure (X2)

        call mma_Allocate(XYZ,3*nAt*8,4)
        iXA0 = 1
        iXA1 = 2
        iXA2 = 3
        iXA3 = 4
        call Expand_Coor(MEP(1,iX0),nAt,XYZ(1,iXA0),mAt)
        call Expand_Coor(MEP(1,iX1),nAt,XYZ(1,iXA1),mAt)
        call Expand_Coor(RP_Centers(1,1,ipX2),nAt,XYZ(1,iXA2),mAt)
        call Expand_Coor(RP_Centers(1,1,ipX3),nAt,XYZ(1,iXA3),mAt)
        if (Invar) then
          call Superpose_w(XYZ(1,iXA0),XYZ(1,iXA2),W,mAt,RMSD,RMax)
          call Fix_Symmetry(XYZ(1,iXA0),nAt,iStab)
          call Superpose_w(XYZ(1,iXA1),XYZ(1,iXA2),W,mAt,RMSD,RMax)
          call Fix_Symmetry(XYZ(1,iXA1),nAt,iStab)
          call Superpose_w(XYZ(1,iXA3),XYZ(1,iXA2),W,mAt,RMSD,RMax)
          call Fix_Symmetry(XYZ(1,iXA3),nAt,iStab)
        end if
        iX0 = iXA0
        iX1 = iXA1
        iX3 = iXA3

        ! Compute deviation=(X1-X0)*(X2-X1).
        ! If it is lower than 0.8, slow down and do linear interpolation
        ! (the direction is probably broken)

        if (iSaddle > 2) then
          call dcopy_(nRP,XYZ(1,iX1),1,Vec(1,1),1)
          call daxpy_(nRP,-One,XYZ(1,iX0),1,Vec(1,1),1)
          call dcopy_(nRP,RP_Centers(1,1,ipX2),1,Vec(1,2),1)
          call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,2),1)
          deviation = dmwdot(nAt,mAt,Vec(1,1),Vec(1,2))
          R11 = dmwdot(nAt,mAt,Vec(1,2),Vec(1,2))
          R22 = dmwdot(nAt,mAt,Vec(1,1),Vec(1,1))
          deviation = deviation/sqrt(R11*R22)
          if (deviation < 0.85d0) then
            quadratic = .false.
            dHSR = SadStep*0.8d0
            Go To 35
          end if
        end if

        ! Compute R11=(X3-X1)**2 and R22=(X2-X1)**2

        call dcopy_(nRP,XYZ(1,iX3),1,Vec(1,1),1)
        call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,1),1)
        R11 = dmwdot(nAt,mAt,Vec(1,1),Vec(1,1))
        call dcopy_(nRP,RP_Centers(1,1,ipX2),1,Vec(1,2),1)
        call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,2),1)
        R22 = dmwdot(nAt,mAt,Vec(1,2),Vec(1,2))
        R1R2 = dmwdot(nAt,mAt,Vec(1,2),Vec(1,1))

        ! The direction of the previous iteration is far from the R-P direction

        tmp = R1R2/(sqrt(R11)*sqrt(R22))
        if (tmp < 0.3d0) then
          quadratic = .false.
          dHSR = SadStep*0.8d0
        elseif (tmp < Zero) then
          quadratic = .false.
          dHSR = SadStep*0.6
        else
          dHSR = SadStep*1.0d0
          if (isaddle > 1) then
            dHSR = SadStep*1.3d0
          end if
        end if
35      continue
        call mma_deallocate(XYZ)
        call mma_deallocate(MEP)
        if (.not. quadratic) call mma_deallocate(Vec)

        ! Slow down close to the TS or in the first iteration

        if (HSR <= 2.5d0*SadStep) then
          dHSR = SadStep*0.55d0
        elseif (quadratic .and. (HSR <= 4.0d0*SadStep)) then
          dHSR = SadStep*0.75d0
        end if
      else
        dHSR = SadStep*0.7d0
      end if
      Delta = dHSR/HSR
      if (Mode == 'P') Delta = One-Delta
      HSR = HSR-dHSR
    end if
    TmpA(6*nAt+5) = Update
    TmpA(6*nAt+4) = HSR
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! This is written here only the first time around

    FindTS = .false.
    HSR0 = HSR

    ! Already findTS!

    FindTS = HSR <= (1.5d0*SadStep)
    if (FindTS) then
      write(6,*) '**************************'
      write(6,*) '* Enable TS optimization *'
      write(6,*) '**************************'
      Update = 2.0d0
      Delta = 0.5d0-6.25d0*(E2-E1)
      Delta = min(Delta,0.75d0)
      Delta = max(Delta,0.25d0)
      if (Mode == 'R') then
        HSR = (One-Delta)*HSR
      else
        HSR = Delta*HSR
      end if

      ! Calculate tangent vector

      call mma_allocate(TanVec,3*nAt,label='TanVec')
      if (Mode == 'R') then
        call Calc_LSTvec(3*nAt,RP_Centers(1,1,2),RP_Centers(1,1,1),TanVec,Invar)
      else
        call Calc_LSTvec(3*nAt,RP_Centers(1,1,1),RP_Centers(1,1,2),TanVec,Invar)
      end if
      call Put_dArray('TanVec',TanVec,3*nAt)
      call mma_deallocate(TanVec)
    else

      ! Do not go too fast the first time

      dHSR = SadStep*0.7d0
      Delta = dHSR/HSR
      HSR = HSR-dHSR
      if (Mode == 'P') Delta = One-Delta
    end if
    call dcopy_(3*nAt,RP_Centers(1,1,1),1,TmpA(1),1)
    call dcopy_(3*nAt,RP_Centers(1,1,2),1,TmpA(1+3*nAt),1)
    TmpA(6*nAt+1) = E1
    TmpA(6*nAt+2) = E2
    TmpA(6*nAt+3) = HSR0
    TmpA(6*nAt+4) = HSR
    TmpA(6*nAt+5) = Update
  end if

  call Put_dArray('Saddle',TmpA,nSaddle)
  call mma_deallocate(TmpA)

  !*********************************************************************
  !                                                                    *
  !                    Some paper work                                 *
  !                                                                    *
  !*********************************************************************
  ! Set the point with the highest energy as the reference
  ! structure. Put the reference structure on the runfile.

  write(6,*)
  write(6,'(A)') ' -- TS optimization a la the Saddle approach'
  if (FindTS) then
    write(6,'(A)') '   Last Macro iteration'
    call Merge_Lists(Mode,nAt)
  end if
  if (Mode == 'R') then
    iRef = 2
    iOpt = 1
    write(6,'(A)') '     Reference structure: product side'
    write(6,'(A,F15.8)') '       Associated Energy: ',E2
    write(6,'(A)') '     Optimized structure: reactant side'
    write(6,'(A,F15.8)') '       Associated Energy: ',E1
  else
    iRef = 1
    iOpt = 2
    write(6,'(A)') '     Reference structure: reactant side'
    write(6,'(A,F15.8)') '       Associated Energy: ',E1
    write(6,'(A)') '     Optimized structure: product side'
    write(6,'(A,F15.8)') '       Associated Energy: ',E2
  end if
  write(6,*)

  ! Align the reference structure with the current structure

  call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
  iRefAlign = 1
  iOptExp = 2
  call Expand_Coor(RP_Centers(1,1,iRef),nAt,XYZ(1,iRefAlign),mAt)
  call Expand_Coor(RP_Centers(1,1,iOpt),nAt,XYZ(1,iOptExp),mAt)
  if (Invar) then
    call Superpose_w(XYZ(1,iRefAlign),XYZ(1,iOptExp),W,mAt,RMSD,RMax)
    call Fix_Symmetry(XYZ(1,iRefAlign),nAt,iStab)
  end if
  call Put_dArray('Ref_Geom',XYZ(1,iRefAlign),3*nAt)
  call mma_deallocate(XYZ)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Fix so that two copies of runfile, vector files, jobiph, etc.
  ! are maintained correctly. This is only done the very first
  ! iteration.

  if (.not. Not_First_Iter) then
    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_Open(LuInput,StdIn)

    ! Make two copies of the virgin runfile

    write(LuInput,'(A)') '> CLONE $Project.RunFile $Project.Reac.RunFile'
    write(LuInput,'(A)') '> CLONE $Project.RunFile $Project.Prod.RunFile'
    write(LuInput,'(A)') '> RM $Project.RunFile'

    ! Make the appropriate links for the other files.

    if (Mode == 'R') then
      write(LuInput,'(A)') '> EXPORT SubProject=.Reac'
    else
      write(LuInput,'(A)') '> EXPORT SubProject=.Prod'
    end if
    write(LuInput,'(A)') '> CLONE $Project.GssOrb $Project$SubProject.GssOrb'
    write(LuInput,'(A)') '> RM -FORCE $Project.GssOrb'

    ! Signal that this is a saddle run, and the first iteration in a branch

    write(LuInput,'(A)') '> EXPORT MOLCAS_SADDLE=1'
    write(LuInput,'(A)') '> EXPORT SADDLE_FIRST=1'
    close(LuInput)
  else
    lRP_Post = .false.
  end if

  !*********************************************************************
  !                                                                    *
  !     Compute the starting structure and put it into the runfile     *
  !                                                                    *
  !*********************************************************************
30 continue
  if (quadratic) then

    ! Quadratic interpolation:
    ! Find the solutions of the system

    R1_2 = R11-Two*R1R2+R22
    Delta = HSR**2*(R11*R22*R1_2+HSR**2*(R1R2**2-R11*R22))
    if (Delta < Zero) then
      write(6,*) 'Delta is negative!!!'
      quadratic = .false.
      Go To 30
    end if
    C = (HSR**2*(Two*R1R2**2-(R11+R1R2)*R22)+(-Two*R1R2+R22)*sqrt(Delta))/(R11*R22*R1_2)
    D = (HSR**2*(R22-R1R2)+sqrt(Delta))/(R22*R1_2)

    ! Compute the next structure

    call daxpy_(nRP,(C-One),Vec(1,1),1,Vec(1,1),1)
    call daxpy_(nRP,D,Vec(1,2),1,Vec(1,1),1)
    call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
    iRA1 = 1
    iRA2 = 2
    call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt)
    call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt)
    if (Mode == 'R') then
      if (Invar) then
        call Superpose_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,RMSD,RMax)
        call Fix_Symmetry(XYZ(1,iRA2),nAt,iStab)
      end if
      call daxpy_(nRP,One,XYZ(1,iRA2),1,Vec(1,1),1)
    else
      if (Invar) then
        call Superpose_w(XYZ(1,iRA1),XYZ(1,iRA2),W,mAt,RMSD,RMax)
        call Fix_Symmetry(XYZ(1,iRA1),nAt,iStab)
      end if
      call daxpy_(nRP,One,XYZ(1,iRA1),1,Vec(1,1),1)
    end if
    call mma_deallocate(XYZ)
    j = 1
    do iCnttp=1,nCnttp
      if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
        do iCnt=1,dbsc(iCnttp)%nCntr
          do i=1,3
            dbsc(iCnttp)%Coor(i,iCnt) = Vec(j,1)
            j = j+1
          end do
        end do
      end if
    end do
    if (Not_First_Iter) call Put_Coord_New(Vec(1,1),nAt)
    call mma_deallocate(Vec)
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Linear interpolation

    call mma_allocate(TmpA,nRP,label='TmpA')
    call DZero(TmpA,nRP)
    call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
    iRA1 = 1
    iRA2 = 2
    call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt)
    call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt)
    if (Invar) then
      if (Mode == 'R') then
        call Superpose_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,RMSD,RMax)
        call Fix_Symmetry(XYZ(1,iRA2),nAt,iStab)
      else
        call Superpose_w(XYZ(1,iRA1),XYZ(1,iRA2),W,mAt,RMSD,RMax)
        call Fix_Symmetry(XYZ(1,iRA1),nAt,iStab)
      end if
    end if
    call daxpy_(nRP,(One-Delta),XYZ(1,iRA1),1,TmpA,1)
    call daxpy_(nRP,(Delta),XYZ(1,iRA2),1,TmpA,1)
    call mma_deallocate(XYZ)
    j = 1
    do iCnttp=1,nCnttp
      if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
        do iCnt=1,dbsc(iCnttp)%nCntr
          do i=1,3
            dbsc(iCnttp)%Coor(i,iCnt) = TmpA(j)
            j = j+1
          end do
        end do
      end if
    end do
    if (Not_First_Iter) call Put_Coord_New(TmpA,nAt)
    call mma_deallocate(TmpA)
  end if

  call mma_deallocate(iStab)
  call mma_deallocate(W)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write constraint on the runfile for slapaf to pick up later.

100 continue
  Lu_UDC = 97
  Lu_UDC = IsFreeUnit(Lu_UDC)
  call Molcas_Open(Lu_UDC,'UDC.Saddle')
  write(Lu_UDC,*) 'R = Sphere'
  write(Lu_UDC,*) 'Value'
  write(Lu_UDC,*) 'R = ',HSR,' soft'
  write(Lu_UDC,*) 'END'
  close(Lu_UDC)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Saddle
