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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               1995, Martin Schuetz                                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NewOrb_SCF(AllowFlip)
!***********************************************************************
!                                                                      *
!     purpose: Diagonalize Fock matrix to get new orbitals             *
!                                                                      *
!     input:                                                           *
!       Fock    : Fock matrix of length nFock                          *
!       CMO     : orthonormal vectors from previous iteration of       *
!                 length nCMO                                          *
!                                                                      *
!     output:                                                          *
!       CMO     : orthonormal vectors in current iteration             *
!       FMOMax  : Max Element of occ/virt block in Fock Matrix,        *
!                 transformed into MO basis                            *
!       EOrb    : orbital energies of length nEOrb                     *
!                                                                      *
!***********************************************************************

use SpinAV, only: Do_SpinAV
use InfSCF, only: MxConstr, DoHLGap, HLGap, FlipThr, FCKAuf, MaxBas, MaxOrf, nnB, WarnCFG, nnFr, Aufb, nSym, TEEE, RotFac, RotLev, &
                  ScrFac, RotMax, MaxBOF, nBas, nBB, nBO, nBT, nConstr, nFro, nOcc, nOrb, TimFld, nD, Iter, Scrmbl, FMOMax
use Constants, only: Zero, One, Two
use stdalloc, only: mma_allocate, mma_deallocate
use SCF_Arrays, only: Ovrlp, EOrb, Fock => FockAO, CMO

implicit none
logical AllowFlip
real*8, dimension(:), allocatable :: eConstr, FckM, FckS, HlfF, TraF, Scratch, Temp
real*8 Fia, GapAdd, EHOMO, ELUMO, Q, Tmp1, Tmp2, Dummy, Tmp, CPU1, CPU2, Tim1, Tim2, Tim3, WhatEver, Tmp0
real*8, external :: DDot_, Random_Molcas
integer, external :: iDAMax_
logical em_On, Scram
integer, dimension(:), allocatable :: iFerm
integer iCMO, iiBT, jEOr, iptr, nOrbmF, nOccmF, nVrt, ia, ij, nsDg, iChk, iiB, iAddGap, iSym, iBas, Ind, iD, iOvlpOff, kConstr, &
        iConstr, nj, iFC, j, jj, ijBas, iOrb, jOrb, IndII, IndJJ, iDum, i, kBas, Muon_I, jBas, IndIJ, jOff, ii, kk, kOff, lConstr, &
        kCMO, keOR, iErr, nFound, Muon_J, iOff, kOrb
integer, save :: iSeed = 13
integer Fermion_Type
logical :: Muons_Present = .false.

!                                                                      *
!***********************************************************************
!                                                                      *
call Timing(Cpu1,Tim1,Tim2,Tim3)
#ifdef _DEBUGPRINT_
call NrmClc(Fock,size(Fock),'NewOrb','Fock')
#endif

Scram = (Scrmbl .and. (iter == 1))
nSdg = 1
if (Do_SpinAV) nSdg = 2
if (MxConstr > 0) call mma_allocate(eConstr,nSdg*MxConstr,Label='eConstr')

if (.not. AllowFlip) then
  if (.not. DoHLgap) then
    DoHLgap = .true.
    HLgap = FlipThr
  end if
end if
! Allocate memory for orbital homeing
if (.not. FckAuf) then
  call mma_allocate(Scratch,MaxBas**2,Label='Scratch')
  call mma_allocate(Temp,MaxBas**2,Label='TempX')
end if
! Allocate memory for modified Fock matrix
call mma_allocate(FckM,nBT,Label='FckM')
! Allocate memory for squared Fock matrix
call mma_allocate(FckS,MaxBas**2,Label='FckSX')
! Allocate memory for half-transformed Fock matrix
call mma_allocate(HlfF,MaxBOF,Label='HlfF')
! Allocate memory for transformed Fock matrix (triangular)
call mma_allocate(TraF,MaxOrF*(MaxOrF+1)/2,Label='TraF')
! Allocate memory for fermi index array
call mma_allocate(iFerm,nBB,Label='iFerm')
call Get_iArray('Fermion IDs',iFerm,nnB)
iChk = 0
do iiB=1,nnB
  iChk = iChk+iFerm(iiB)
end do
em_On = ((iChk /= 0) .and. (iChk /= nnB))

FMOMax = Zero
WarnCfg = .false.
do iD=1,nD

  ! Modify Fock matrix
  call dcopy_(nBT,Fock(1,iD),1,FckM,1)
  if (nnFr > 0) call ModFck(FckM,Ovrlp,nBT,CMO(1,iD),nBO,nOcc(1,iD))
  ! Prediagonalize Fock matrix
  iAddGap = 0
  GapAdd = Zero
  if ((.not. Aufb) .and. DoHLgap) then
    ij = 1
    iCMO = 1
    jEOr = 1
    Ehomo = -1.0d6
    Elumo = 1.0d6
    do iSym=1,nSym
      iiBT = nBas(iSym)*(nBas(iSym)+1)/2
      nOrbmF = nOrb(iSym)-nFro(iSym)
      nOccmF = nOcc(iSym,iD)-nFro(iSym)
      nVrt = nOrb(iSym)-nOcc(iSym,iD)
      iCMO = iCMO+nBas(iSym)*nFro(iSym)
      jEOr = jEOr+nFro(iSym)
      call Square(FckM(ij),FckS,1,nBas(iSym),nBas(iSym))
      if (nOccmF > 0) then
        call DGEMM_('N','N',nBas(iSym),nOccmF,nBas(iSym), &
                    One,FckS,nBas(iSym), &
                    CMO(iCMO,iD),nBas(iSym), &
                    Zero,HlfF,nBas(iSym))
        call DGEMM_Tri('T','N',nOccmF,nOccmF,nBas(iSym), &
                       One,CMO(iCMO,iD),nBas(iSym), &
                       HlfF,nBas(iSym), &
                       Zero,TraF,nOccmF)
#       ifdef _DEBUGPRINT_
        call Triprt('Occupied Fock matrix in MO basis','(20F10.4)',TraF,nOccmF)
#       endif
        nOccmF = nOccmF-nConstr(iSym)
        call NIdiag(TraF,CMO(iCMO,iD),nOccmF,nBas(iSym))
        nOccmF = nOccmF+nConstr(iSym)
#       ifdef _DEBUGPRINT_
        call Triprt('Occupied Fock matrix in MO basis','(20F10.4)',TraF,nOccmF)
#       endif
        do iBas=1,nOccmF
          ind = iBas*(iBas+1)/2
          Ehomo = max(Ehomo,TraF(ind))
        end do
      end if

      iCMO = iCMO+nOccmF*nBas(iSym)
      jEOr = jEOr+nOccmF
      if (Do_SpinAV) then
        nVrt = nVrt-nConstr(iSym)
        iCMO = iCMO+nConstr(iSym)*nBas(iSym)
        jEOr = jEOr+nConstr(iSym)
      end if
      if (nVrt > 0) then
        call DGEMM_('N','N',nBas(iSym),nVrt,nBas(iSym), &
                    One,FckS,nBas(iSym), &
                    CMO(iCMO,iD),nBas(iSym), &
                    Zero,HlfF,nBas(iSym))
        call DGEMM_Tri('T','N',nVrt,nVrt,nBas(iSym), &
                       One,CMO(iCMO,iD),nBas(iSym), &
                       HlfF,nBas(iSym), &
                       Zero,TraF,nVrt)
#       ifdef _DEBUGPRINT_
        call Triprt('Virtual Fock matrix in MO basis','(20F10.4)',TraF,nVrt)
#       endif
        call NIdiag(TraF,CMO(iCMO,iD),nVrt,nBas(iSym))
#       ifdef _DEBUGPRINT_
        call Triprt('Virtual Fock matrix in MO basis','(20F10.4)',TraF,nVrt)
#       endif
        do iBas=1,nVrt
          ind = iBas*(iBas+1)/2
          Elumo = min(Elumo,TraF(ind))
        end do
      end if
      if (Do_SpinAV) then
        nVrt = nVrt+nConstr(iSym)
        iCMO = iCMO-nConstr(iSym)*nBas(iSym)
        jEOr = jEOr-nConstr(iSym)
      end if
      iCMO = iCMO+nVrt*nBas(iSym)
      jEOr = jEOr+nVrt
      ij = ij+iiBT
    end do
#   ifdef _DEBUGPRINT_
    write(6,'(a,F12.6)') 'E(homo)   ',Ehomo
    write(6,'(a,F12.6)') 'E(lumo)   ',Elumo
    write(6,'(a,F12.6)') 'E(gap)    ',Elumo-Ehomo
#   endif
    WarnCfg = ((Elumo-Ehomo < Zero) .or. WarnCfg)
    if (Elumo-Ehomo < HLgap) then
      iAddGap = 1
      GapAdd = HLgap-Elumo+Ehomo
#     ifdef _DEBUGPRINT_
      write(6,'(a,F12.6)') 'E(add)    ',GapAdd
#     endif
    end if
  end if
  ! Diagonalize Fock matrix in non-frozen molecular basis
  ij = 1
  iCMO = 1
  jEOr = 1
  iOvlpOff = 1
  do iSym=1,nSym
    iiBT = nBas(iSym)*(nBas(iSym)+1)/2
    nOrbmF = nOrb(iSym)-nFro(iSym)
    nOccmF = nOcc(iSym,iD)-nFro(iSym)
    nVrt = nOrb(iSym)-nOcc(iSym,iD)
    ! Find the proper pointers to CMO and EOr
    iCMO = iCMO+nBas(iSym)*nFro(iSym)
    jEOr = jEOr+nFro(iSym)
    if (nOrbmF > 0) then
      call Square(FckM(ij),FckS,1,nBas(iSym),nBas(iSym))
      ! Transform Fock matrix to the basis from previous iteration
      call DGEMM_('N','N',nBas(iSym),nOrbmF,nBas(iSym), &
                  One,FckS,nBas(iSym), &
                  CMO(iCMO,iD),nBas(iSym), &
                  Zero,HlfF,nBas(iSym))
      call DGEMM_Tri('T','N',nOrbmF,nOrbmF,nBas(iSym), &
                     One,CMO(iCMO,iD),nBas(iSym), &
                     HlfF,nBas(iSym), &
                     Zero,TraF,nOrbmF)

      ! Constrained SCF section begins
      kConstr = 1
      do iConstr=nConstr(iSym),1,-1
        nj = nOccmF-iConstr
        ifc = 1+nj*(nj+1)/2
        eConstr(kConstr) = TraF(ifc+nj)
        call FZero(TraF(ifc),nj)
        do j=nj+1,nOrbmF-1
          jj = 1+j*(j+1)/2+nj
          Traf(jj) = Zero
        end do
        Traf(ifc+nj) = -0.666d6*dble(iConstr) ! for sorting
        kConstr = kConstr+1
      end do
      if (Do_SpinAV) then
        do iConstr=0,nConstr(iSym)-1
          nj = nOccmF+iConstr
          ifc = 1+nj*(nj+1)/2
          eConstr(kConstr) = TraF(ifc+nj)
          call FZero(TraF(ifc),nj)
          do j=nj+1,nOrbmF-1
            jj = 1+j*(j+1)/2+nj
            TraF(jj) = Zero
          end do
          TraF(ifc+nj) = 0.666d6*dble(kConstr) ! for sorting
          kConstr = kConstr+1
        end do
      end if
      ! Constrained SCF section ends

      ! get max element of Fock matrix
      !do iBas=2,nOrbmF
      !  do jBas=1,iBas-1
      !    ijBas = iBas*(iBas-1)/2+jBas-1+1
      !    FMOMax = Max(Abs(TraF(ijBas)),FMOMax)
      !  end do
      !end do
      ! get max element of occ/virt Block of Fock matrix

      if (Teee) then
        do iBas=2,nOrbmF
          do jBas=1,iBas-1
            ijBas = iBas*(iBas-1)/2+jBas
            FMOMax = max(abs(TraF(ijBas)),FMOMax)
          end do
        end do
      else if ((nOccmF > 0) .and. (nVrt > 0)) then
        iptr = 1+nOccmF*(nOccmF+1)/2
        do ia=1,nVrt
          Fia = abs(TraF(iptr+IDAMAX_(nOccmF,TraF(iptr),1)-1))
          FMOMax = max(Fia,FMOMax)
          iptr = iptr+nOccmF+ia
        end do
      end if

      ! Modify Fock matrix to enhance convergence

      ! Add to homo lumo gap

      if (iAddGap == 1) then
        do iOrb=nOccmF+1,nOrbmF
          ind = iOrb*(iOrb+1)/2
          TraF(ind) = TraF(ind)+GapAdd
        end do
      end if
      ind = 1
      do iOrb=1,nOrbmF
        do jOrb=1,iOrb
          ! Scale OV elements of fock matrix
          if ((iOrb > nOcc(iSym,iD)) .and. (jOrb <= nOcc(iSym,iD))) TraF(ind) = RotFac*TraF(ind)
          ! Levelshift virtual diagonal matrix elements
          if ((iOrb > nOcc(iSym,iD)) .and. (iOrb == jOrb)) TraF(ind) = TraF(ind)+RotLev
          ind = ind+1
        end do
      end do
      ! Add scrambling
      if (Scram) then
        ind = 1
        do iOrb=1,nOrbmF
          do jOrb=1,iOrb
            if (iOrb /= jOrb) then
              q = ScrFac*(Two*Random_Molcas(iSeed)-One)
              TraF(ind) = TraF(ind)+q
            end if
            ind = ind+1
          end do
        end do
      end if

      ! Scale OV elements if too big.
      indii = 1
      indij = 1
      do iOrb=1,nOrbmF
        indjj = 1
        do jOrb=1,iOrb
          if ((iOrb > nOcc(iSym,iD)) .and. (jOrb <= nOcc(iSym,iD))) then
            tmp1 = max(abs(TraF(indii)-TraF(indjj)),1.0d-3)
            tmp2 = abs(TraF(indij)/tmp1)
            if ((tmp2 > RotMax) .and. (abs(TraF(indij)) > 0.001d0)) TraF(indij) = TraF(indij)*RotMax/tmp2
          end if
          indij = indij+1
          indjj = indjj+jOrb+1
        end do
        indii = indii+iOrb+1
      end do

      ! Constrained SCF section begins
      kConstr = 1
      do iConstr=nConstr(iSym),1,-1
        nj = nOccmF-iConstr
        ifc = 1+nj*(nj+1)/2
        call FZero(TraF(ifc),nj)
        do j=nj+1,nOrbmF-1
          jj = 1+j*(j+1)/2+nj
          TraF(jj) = Zero
        end do
        TraF(ifc+nj) = -0.666d6*dble(iConstr) ! for sorting
        kConstr = kConstr+1
      end do
      if (Do_SpinAV) then
        do iConstr=0,nConstr(iSym)-1
          nj = nOccmF+iConstr
          ifc = 1+nj*(nj+1)/2
          call FZero(TraF(ifc),nj)
          do j=nj+1,nOrbmF-1
            jj = 1+j*(j+1)/2+nj
            TraF(jj) = Zero
          end do
          TraF(ifc+nj) = 0.666d6*dble(kConstr) ! for sorting
          kConstr = kConstr+1
        end do
      end if
      ! Constrained SCF section ends

#     ifdef _DEBUGPRINT_
      call Triprt('Fock matrix in MO basis after modification','(10F10.4)',TraF,nOrbmF)
#     endif

      ! Diagonalize and form orbital energies

      Dummy = Zero
      iDum = 0

      ! Store the original CMOs for root following.

      call DCopy_(nBas(iSym)**2,CMO(iCMO,iD),1,FckS,1)

      call Diag_Driver('V','A','L',nOrbmF,TraF,TraF,nOrbmF,Dummy,Dummy,iDum,iDum,EOrb(jEOr,iD),CMO(iCMO,iD),nBas(iSym),0,-1,'J', &
                       nFound,iErr)

      ! Fix standard phase of the orbitals

      do i=1,nBas(iSym)
        call VecPhase(CMO(iCMO+(i-1)*nBas(iSym),iD),nBas(iSym))
      end do
#     ifdef _DEBUGPRINT_
      call NrmClc(Fcks,nbas(iSym)*nOrb(iSym),'NewOrb','Old CMOs')
      call NrmClc(CMO(iCMO,iD),nbas(iSym)*nOrb(iSym),'NewOrb','New CMOs')
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Reorder the orbitals to preserve e/m partition.

      if (.not. em_On) Go To 110 ! skip if all orbitals are of
      ! the same type.

      do iOrb=1,nOrb(iSym)-1    ! Loop over the old orbitals

        ! Compute a check sum which is not zero if the orbital
        ! is a muonic orbital.

        tmp = Zero
        do kBas=0,nBas(iSym)-1
          tmp = tmp+dble(iFerm(jEOr+kBas))*abs(FckS((iOrb-1)*nBas(iSym)+kBas+1))
        end do
        Muon_i = 0                  ! electronic
        if (tmp /= Zero) Muon_i = 1! muonic
        Muons_Present = (Muons_Present .or. (Muon_i == 1))
        !write(6,*) 'iOrb,Muon_i,tmp=',iOrb,Muon_i,tmp

        ! Loop over the new orbitals and test if it is of the
        ! same type. i.e. fermionic or electronic.

        do jOrb=iOrb,nOrb(iSym)
          tmp = Zero
          do kBas=0,nBas(iSym)-1
            tmp = tmp+dble(iFerm(jEOr+kBas))*abs(CMO(iCMO+(jOrb-1)*nBas(iSym)+kBas,iD))
          end do

          Muon_j = 0                  ! electronic
          if (tmp /= Zero) Muon_j = 1 ! muonic

          ! If orbital and fermion index are identical fine.

          if ((iOrb == jOrb) .and. (Muon_i == Muon_j)) Go To 678

          ! If fermion index the same swap orbital and the
          ! corresponding orbital energy in the new list.

          if (Muon_i == Muon_j) then
            Tmp = EOrb(jEOr-1+iOrb,iD)
            EOrb(jEOr-1+iOrb,iD) = EOrb(jEOr-1+jOrb,iD)
            EOrb(jEOr-1+jOrb,iD) = Tmp
            call DSwap_(nBas(iSym),CMO(iCMO+(iOrb-1)*nBas(iSym),iD),1,CMO(iCMO+(jOrb-1)*nBas(iSym),iD),1)
            Go To 678
          end if

        end do   !  jOrb

        ! Arrive at this point when the iOrb'th orbital in the
        ! new list is of the same type as in the old list.

678     continue

      end do      !  iOrb
110   continue
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Order the orbitals in the same way as previously, that is,
      ! do not populate according to the aufbau principle.

      if (.not. FckAuf) then
        !write(6,*) 'Follow the orbitals'

        ! Form  C^+ S  C for the old orbitals

        call FZero(Scratch,nOrb(iSym)*nBas(iSym))
        call Square(Ovrlp(iOvlpOff),Temp,1,nBas(iSym),nBas(iSym))
        call DGEMM_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym), &
                    One,FckS,nBas(iSym), &
                    Temp,nBas(iSym), &
                    Zero,Scratch,nOrb(iSym))

        Fermion_type = 1
        if (Muons_present) Fermion_type = 0
        do iOrb=1,nOrb(iSym)-1  ! Loop over the old orbitals

          ! Don't apply non-aufbau principle for the electrons!

          iOff = (iOrb-1)*nBas(iSym)+iCMO

          ! Identify orbital type.

          tmp = Zero
          do kBas=0,nBas(iSym)-1
            tmp = tmp+dble(iFerm(jEOr+kBas))*abs(FckS((iOrb-1)*nBas(iSym)+kBas+1))
          end do

          Muon_i = 0                  ! electronic
          if (tmp /= Zero) Muon_i = 1 ! muonic
          if (Muon_i == Fermion_Type) cycle

          !write(6,*) 'iOrb,Muon_i=',iOrb,Muon_i

          kOrb = 0
          Tmp0 = Zero
          do jOrb=1,nOrb(iSym)   ! Loop over the new orbitals
            jOff = (jOrb-1)*nBas(iSym)+iCMO

            tmp = Zero
            do kBas=0,nBas(iSym)-1
              tmp = tmp+dble(iFerm(jEOr+kBas))*abs(FckS((jOrb-1)*nBas(iSym)+kBas+1))
            end do

            Muon_j = 0                  ! electronic
            if (tmp /= Zero) Muon_j = 1 ! muonic
            if (Muon_j /= Muon_i) cycle

            Tmp1 = abs(DDot_(nBas(iSym),Scratch(iOrb),nOrb(iSym),CMO(jOff,iD),1))
            if (Tmp1 > Tmp0) then
              Tmp0 = Tmp1
              kOrb = jOrb
            end if
          end do

          !write(6,*) 'Fermion_Type=',Fermion_type
          !write(6,*) 'kOrb,Tmp0=',kOrb,Tmp0
          !write(6,*) 'Muon_i, Muon_j=',Muon_i, Muon_j

          if (iOrb /= kOrb) then
            ii = iOrb+jEOr-1
            kk = kOrb+jEOr-1
            tmp = EOrb(ii,iD)
            EOrb(ii,iD) = EOrb(kk,iD)
            EOrb(kk,iD) = tmp
            kOff = (kOrb-1)*nBas(iSym)+iCMO
            call DSwap_(nBas(iSym),CMO(iOff,iD),1,CMO(kOff,iD),1)
          end if
        end do
#       ifdef _DEBUGPRINT_
        call NrmClc(CMO(iCMO,iD),nbas(iSym)*nOrb(iSym),'NewOrb','New CMOs')
#       endif
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *

      ! Constrained SCF section begins

      if (nConstr(iSym) > 0) then
        do kConstr=1,nConstr(iSym)
          iConstr = jEOr-1+kConstr
          EOrb(iConstr,iD) = 0.666d6*dble(kConstr)
        end do
        call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),nOccmF,nBas(iSym),1,.true.)
        iConstr = 1
        do kConstr=nConstr(iSym),1,-1
          lConstr = jEOr+nOccmF-kConstr
          EOrb(lConstr,iD) = eConstr(iConstr)
          iConstr = iConstr+1
        end do
        if (Do_SpinAV) then
          do kConstr=nConstr(iSym),1,-1
            iConstr = jEOr+nOrbmF-kConstr
            EOrb(iConstr,iD) = -0.666d6*dble(kConstr)
          end do
          kCMO = iCMO+nOccmF*nBas(iSym)
          kEOr = jEOr+nOccmF
          call SortEig(EOrb(kEOr,iD),CMO(kCMO,iD),nVrt,nBas(iSym),1,.true.)
          iConstr = nConstr(iSym)+1
          do kConstr=1,nConstr(iSym)
            lConstr = kEOr+kConstr
            EOrb(lConstr,iD) = eConstr(iConstr)
            iConstr = iConstr+1
          end do
        end if
      end if
      ! Constrained SCF section ends
    end if
    ! Update pointers
    iCMO = iCMO+nOrbmF*nBas(iSym)
    jEOr = jEOr+nOrbmF
    ij = ij+iiBT
    iOvlpOff = iOvlpOff+nBas(iSym)*(nBas(iSym)+1)/2
  end do

  ! Check orthogonality
  call ChkOrt(iD,Whatever)

end do

! Deallocate memory
call mma_deallocate(iFerm)
call mma_deallocate(TraF)
call mma_deallocate(HlfF)
call mma_deallocate(FckS)
call mma_deallocate(FckM)
if (.not. FckAuf) then
  call mma_deallocate(Scratch)
  call mma_deallocate(Temp)
end if
if (MxConstr > 0) call mma_deallocate(eConstr)

#ifdef _DEBUGPRINT_
do iD=1,nD
  iOff = 1
  jOff = 1
  do iSym=1,nSym
    call RecPrt('CMO',' ',CMO(jOff,iD),nBas(iSym),nOrb(iSym))
    iOff = iOff+nOrb(iSym)
    jOff = jOff+nBas(iSym)*nOrb(iSym)
  end do
end do
call RecPrt('NewOrb_scf: EOr',' ',EOrb,1,size(EOrb))
#endif

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(8) = TimFld(8)+(Cpu2-Cpu1)

return

end subroutine NewOrb_SCF
