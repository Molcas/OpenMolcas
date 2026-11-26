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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> How far the active space moved between the starting orbitals and the
!> converged ones. For every irrep, M = C_in^T S C_out (C_in: active
!> orbitals the macro iterations start from, C_out: converged active
!> orbitals, S: AO overlap). Sum_ij M(i,j)^2 = Tr(P_in P_out) is the part
!> of the start active space that the converged one still spans, and the
!> deficit D = n_act - Sum_ij M(i,j)^2 counts the orbitals' worth of the
!> start that left the active space. D and the cosines of the principal
!> angles (the singular values of M) do not depend on rotations within the
!> active spaces; both orbital sets must be expanded in the same basis
!> (in a basis-set series, the start projected by MORT). The weight
!> outside the converged active space of start orbital i,
!> w_i = 1-Sum_j M(i,j)^2, tells which orbital left (Sum_i w_i = D).
!> The active space counts as lost when the total deficit over all irreps
!> exceeds ThrAct (keyword ACTThreshold, default 1.0; 0 switches the check
!> off); RASSCF then ends with _RC_GENERAL_ERROR_, unless EXPERT is given.
!> The check is done only in a restart, when the starting orbitals hold the
!> wave function of this calculation (ReadVc, Act_Set_Restart, Act_Same_WF:
!> the same active electrons, spin, symmetry, active orbitals per irrep,
!> roots and weights): a JOBIPH/JOBOLD or RASSCF HDF5 file of the same wave
!> function, the RASSCF orbitals of the runfile after the first iteration
!> of a geometry optimization or molecular dynamics, and an HDF5 file
!> projected by MORT (basis-set series) with the same active orbitals per
!> irrep. Not from SCF orbitals, computed guesses or text orbital files
!> (no wave-function information), nor when the orbitals of another state
!> or active space are reused: then the CASSCF is expected to change the
!> active orbitals.
!> A second check (warning only, also without a restart) looks at the
!> occupation numbers of the converged active natural orbitals: an orbital
!> with occupation above 2-ThrOcc or below ThrOcc (keyword OCCThreshold,
!> default 0.005; 0 switches it off) adds hardly any correlation (e.g. a
!> lone pair that left the correlated space as a pure atomic orbital). The
!> flagged orbitals are printed with their largest atomic contributions
!> (Mulliken gross populations). If the starting orbitals carry occupation
!> numbers (JOBIPH/JOBOLD, HDF5, orbital file), the start and final
!> occupations are printed side by side: for each converged natural orbital
!> its occupation in the starting active density, with the starting active
!> orbital it overlaps most.
module Active_Space_Check

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
private

! the threshold on the total deficit (ACTThreshold) and whether a lost
! active space stops the calculation (not with EXPERT)
real(kind=wp) :: ThrAct = 1.0_wp
logical(kind=iwp) :: ActStop = .true.
! the threshold of the occupation check (OCCThreshold)
real(kind=wp) :: ThrOcc = 0.005_wp
! copy of the starting orbitals (Act_Check_Start), the check applies
! (kept for the cosines in a restart, or to match the occupations of the
! start)
real(kind=wp), allocatable :: CMO_Start(:), Occ_Start(:)
logical(kind=iwp) :: Have_Start = .false.
! the cosines are computed (restart with the same active space)
logical(kind=iwp) :: Do_Cos = .false.
! the starting orbitals carry occupation numbers
logical(kind=iwp) :: Have_Occ_Start = .false.
! the starting orbitals hold the active space of an earlier CASSCF/RASSCF
logical(kind=iwp) :: Restart = .false.
! active orbitals per irrep of the start file (type indices; -1: unknown)
integer(kind=iwp) :: nAct_File(8) = -1

public :: Act_Check_Final, Act_Check_Free, Act_Check_Start, Act_In_Loop, Act_Same_WF, Act_Set_Restart, Act_Warn_SCF, ActStop, &
          ThrAct, ThrOcc

contains

!> One line of warning when the starting orbitals come from SCF (runfile,
!> SCF HDF5 file, or an orbital file with the SCF title): their active
!> orbitals are only those at the Fermi level of the SCF.
subroutine Act_Warn_SCF()
  write(u6,'(6X,A)') 'WARNING: RASSCF starts from SCF orbitals; check that the active orbitals are the intended ones.'
end subroutine Act_Warn_SCF

!> Whether the starting orbitals hold the active space of an earlier
!> CASSCF/RASSCF (set by ReadVc); TInd: the type indices of the file
!> (tpstr2tpidx: 3, 4, 5 = RAS1, RAS2, RAS3), if it has them. A file with
!> another number of active orbitals in some irrep is not a restart of
!> this active space.
subroutine Act_Set_Restart(L,nTInd,TInd)
  use Orbital_Spaces, only: nBas, nSym
  logical(kind=iwp), intent(in) :: L
  integer(kind=iwp), intent(in), optional :: nTInd, TInd(*)
  integer(kind=iwp) :: iOff, iS, k
  Restart = L
  nAct_File(:) = -1
  if (present(nTInd) .and. present(TInd)) then
    if (nTInd /= sum(nBas(1:nSym))) return
    ! (no type indices on the file)
    if (all(TInd(1:nTInd) <= 0)) return
    iOff = 0
    do iS=1,nSym
      nAct_File(iS) = 0
      do k=iOff+1,iOff+nBas(iS)
        if ((TInd(k) >= 3) .and. (TInd(k) <= 5)) nAct_File(iS) = nAct_File(iS)+1
      end do
      iOff = iOff+nBas(iS)
    end do
  end if
end subroutine Act_Set_Restart

!> The wave function of the start file is the one of this calculation (a
!> restart): the same number of active electrons, spin, state symmetry,
!> active orbitals per irrep, and the same roots with the same weights.
function Act_Same_WF(nActEl_f,iSpin_f,LSym_f,nAsh_f,nRoots_f,iRoot_f,Weight_f) result(Same)
  use CI_Input, only: iSpin, nActEl, STSym
  use Orbital_Spaces, only: nAsh, nSym
  use rasscf_global, only: iRoot, nRoots, Weight
  logical(kind=iwp) :: Same
  integer(kind=iwp), intent(in) :: nActEl_f, iSpin_f, LSym_f, nAsh_f(*), nRoots_f, iRoot_f(*)
  real(kind=wp), intent(in) :: Weight_f(*)
  Same = (nActEl_f == nActEl) .and. (iSpin_f == iSpin) .and. (LSym_f == STSym) .and. (nRoots_f == nRoots)
  if (.not. Same) return
  Same = all(nAsh_f(1:nSym) == nAsh(1:nSym)) .and. all(iRoot_f(1:nRoots) == iRoot(1:nRoots)) .and. &
         all(abs(Weight_f(1:nRoots)-Weight(1:nRoots)) < 1.0e-10_wp)
end function Act_Same_WF

!> A later step of a loop of the input (DO WHILE: geometry optimization,
!> molecular dynamics), MOLCAS_ITER > 1: every step runs the same input,
!> so the orbitals of the step before hold the same wave function.
function Act_In_Loop() result(In_Loop)
  logical(kind=iwp) :: In_Loop
  integer(kind=iwp) :: Iter, iStat
  character(len=80) :: EnvStr
  In_Loop = .false.
  call get_environment_variable('MOLCAS_ITER',EnvStr)
  read(EnvStr,*,iostat=iStat) Iter
  if (iStat == 0) In_Loop = (Iter > 1)
end function Act_In_Loop

!> Keep the orbitals the macro iterations start from: in a restart (the
!> cosines), or when they carry occupation numbers (Occ: as read by ReadVc,
!> zero if the source has none).
subroutine Act_Check_Start(nCMO,CMO,nOcc,Occ)

  use Orbital_Spaces, only: nAsh, nSym
  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nCMO, nOcc
  real(kind=wp), intent(in) :: CMO(nCMO), Occ(nOcc)

  call Free_Start()
  Do_Cos = (ThrAct > Zero) .and. Restart
  ! a start file with another active space is not a restart of this one
  if (Do_Cos .and. (nAct_File(1) >= 0)) then
    if (any(nAct_File(1:nSym) /= nAsh(1:nSym))) then
      write(u6,'(6X,A)') 'Active-space check skipped: the starting orbitals have another number of active orbitals'
      Do_Cos = .false.
    end if
  end if
  Have_Occ_Start = (ThrOcc > Zero) .and. any(Occ(:) /= Zero)
  if (.not. (Do_Cos .or. Have_Occ_Start)) return
  call mma_allocate(CMO_Start,nCMO,Label='CMO_Start')
  CMO_Start(:) = CMO(:)
  if (Have_Occ_Start) then
    call mma_allocate(Occ_Start,nOcc,Label='Occ_Start')
    Occ_Start(:) = Occ(:)
  end if
  Have_Start = .true.

end subroutine Act_Check_Start

!> At the end of the calculation: in a restart, the deficit and the cosines
!> between the converged active space and the starting one (Lost if the
!> total deficit exceeds ThrAct); always (ThrOcc > 0), the occupation check of the
!> converged active natural orbitals. Both use the average natural orbitals
!> and their occupation numbers on JOBIPH (OUTCTL leaves those of the
!> single roots in the arrays of RASSCF). The results are stored on the
!> HDF5 file, if open.
subroutine Act_Check_Final(nCMO,nOcc,Have_H5,H5Id,Lost)

  use OneDat, only: sNoNuc, sNoOri
  use Index_Functions, only: nTri_Elem
  use Orbital_Spaces, only: nAsh, nBas, nFro, nIsh, nSym, nTot1
  use rasscf_files, only: JOBIPH
  use rasscf_global, only: BName, IADR15
  use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _HDF5_
  use mh5, only: mh5_init_attr
#endif

  integer(kind=iwp), intent(in) :: nCMO, nOcc, H5Id
  logical(kind=iwp), intent(in) :: Have_H5
  logical(kind=iwp), intent(out) :: Lost
  integer(kind=iwp) :: i, iA, iAd, iC, iComp, iO, iOpt, iRc, iS, iSyLbl, j, k, nA, nAct, nB, nFlag, nOff1, nOffB, nOffC
  real(kind=wp) :: Def, Def_Irr(8), W
  logical(kind=iwp) :: Do_Match, Have_S
  integer(kind=iwp), allocatable :: Flags(:), iMatch(:), iSym_A(:), iOrb_A(:)
  real(kind=wp), allocatable :: CIn(:,:), CMO(:), COut(:,:), Cos_All(:), M(:,:), Occ(:), Occ_Fin(:), Occ_Match(:), &
                                Ovl(:), S1(:), SC(:,:), Sq(:,:), SV(:)
  character(len=8) :: Label
  character(len=48), allocatable :: Contrib(:)
#include "macros.fh"

  Lost = .false.
#ifndef _HDF5_
  unused_var(H5Id)
  unused_var(Have_H5)
#endif
  nAct = sum(nAsh(1:nSym))
  if ((nAct == 0) .or. ((.not. Do_Cos) .and. (ThrOcc <= Zero))) then
    call Free_Start()
    return
  end if

  ! the AO overlap (symmetry blocked, triangular)
  call mma_allocate(S1,nTot1,Label='S1')
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSyLbl = 1
  Label = 'Mltpl  0'
  call RdOne(iRc,iOpt,Label,iComp,S1,iSyLbl)
  Have_S = (iRc == 0)
  if ((.not. Have_S) .and. Do_Cos) then
    write(u6,'(/,6X,A)') 'Active-space check skipped: no AO overlap on ONEINT'
    Do_Cos = .false.
  end if
  Do_Match = Have_S .and. Have_Occ_Start

  ! the average natural orbitals and their occupation numbers
  call mma_allocate(CMO,nCMO,Label='CMO')
  call mma_allocate(Occ,nOcc,Label='Occ')
  iAd = IADR15(2)
  call DDaFile(JOBIPH,2,CMO,nCMO,iAd)
  call DDaFile(JOBIPH,2,Occ,nOcc,iAd)

  call mma_allocate(Cos_All,nAct,Label='Cos_All')
  call mma_allocate(Flags,nAct,Label='Flags')
  call mma_allocate(iSym_A,nAct,Label='iSym_A')
  call mma_allocate(iOrb_A,nAct,Label='iOrb_A')
  call mma_allocate(iMatch,nAct,Label='iMatch')
  call mma_allocate(Occ_Fin,nAct,Label='Occ_Fin')
  call mma_allocate(Occ_Match,nAct,Label='Occ_Match')
  call mma_allocate(Ovl,nAct,Label='Ovl')
  call mma_allocate(Contrib,nAct,Label='Contrib')
  Flags(:) = 0
  iMatch(:) = 0
  Occ_Match(:) = Zero
  Ovl(:) = Zero
  Contrib(:) = ''
  if (Do_Cos) then
    write(u6,*)
    write(u6,'(6X,A)') 'Active space: overlap of the converged active orbitals with the starting ones'
    write(u6,'(6X,A)') '(deficit = orbitals'' worth of the start that left the active space; 0 = unchanged)'
  end if
  nOff1 = 0
  nOffB = 0
  nOffC = 0
  iA = 0
  iC = 0
  Def_Irr(:) = Zero
  do iS=1,nSym
    nB = nBas(iS)
    nA = nAsh(iS)
    if (nA > 0) then
      iO = nOffC+(nFro(iS)+nIsh(iS))*nB
      call mma_allocate(COut,nB,nA,Label='COut')
      do j=1,nA
        COut(:,j) = CMO(iO+(j-1)*nB+1:iO+j*nB)
      end do
      if (Have_S) then
        call mma_allocate(Sq,nB,nB,Label='Sq')
        call Square(S1(nOff1+1),Sq,1,nB,nB)
        call mma_allocate(SC,nB,nA,Label='SC')
        call DGEMM_('N','N',nB,nA,nB,One,Sq,nB,COut,nB,Zero,SC,nB)
      end if
      if (Have_S .and. Have_Start) then
        call mma_allocate(CIn,nB,nA,Label='CIn')
        do j=1,nA
          CIn(:,j) = CMO_Start(iO+(j-1)*nB+1:iO+j*nB)
        end do
        call mma_allocate(M,nA,nA,Label='M')
        call DGEMM_('T','N',nA,nA,nB,One,CIn,nB,SC,nB,Zero,M,nA)
        call mma_deallocate(CIn)
      end if
      if (Do_Cos) then
        ! weights outside the converged active space
        write(u6,'(/,6X,A,I2,A,I4,A)') 'Irrep',iS,':',nA,' active orbitals'
        write(u6,'(6X,A)') 'start orbital   weight outside the converged active space'
        do i=1,nA
          W = One-sum(M(i,:)**2)
          Def_Irr(iS) = Def_Irr(iS)+W
          if (W > Half) then
            write(u6,'(6X,I8,6X,F10.4,A)') nFro(iS)+nIsh(iS)+i,W,'  <-- left the active space'
          else
            write(u6,'(6X,I8,6X,F10.4)') nFro(iS)+nIsh(iS)+i,W
          end if
        end do
        write(u6,'(6X,A,F9.4)') 'deficit',Def_Irr(iS)
        call mma_allocate(SV,nA,Label='SV')
        call SVD_Values(nA,M,SV)
        write(u6,'(6X,A)') 'cosines of the principal angles:'
        write(u6,'(6X,8F9.4)') (SV(k),k=1,nA)
        Cos_All(iC+1:iC+nA) = SV(:)
        iC = iC+nA
        call mma_deallocate(SV)
      end if
      ! occupation numbers of the converged active natural orbitals
      do j=1,nA
        k = nOffB+nFro(iS)+nIsh(iS)+j
        iSym_A(iA+j) = iS
        iOrb_A(iA+j) = nFro(iS)+nIsh(iS)+j
        Occ_Fin(iA+j) = Occ(k)
        if (ThrOcc > Zero) then
          if (Occ(k) > Two-ThrOcc) Flags(iA+j) = 1
          if (Occ(k) < ThrOcc) Flags(iA+j) = -1
        end if
        if ((Flags(iA+j) /= 0) .and. Have_S) &
          call Atom_Contrib(nB,COut(:,j),SC(:,j),BName(nOffB+1:nOffB+nB),Contrib(iA+j))
        if (Do_Match) then
          ! the occupation of this natural orbital in the starting active
          ! density, sum_i M(i,j)^2 n_i / sum_i M(i,j)^2, and the starting
          ! active orbital it overlaps most
          W = sum(M(:,j)**2)
          if (W > 1.0e-8_wp) Occ_Match(iA+j) = sum(M(:,j)**2*Occ_Start(nOffB+nFro(iS)+nIsh(iS)+1:nOffB+nFro(iS)+nIsh(iS)+nA))/W
          i = maxloc(M(:,j)**2,1)
          iMatch(iA+j) = nFro(iS)+nIsh(iS)+i
          Ovl(iA+j) = M(i,j)**2
        end if
      end do
      iA = iA+nA
      call mma_deallocate(M,safe='*')
      call mma_deallocate(SC,safe='*')
      call mma_deallocate(Sq,safe='*')
      call mma_deallocate(COut)
    end if
    nOff1 = nOff1+nTri_Elem(nB)
    nOffB = nOffB+nB
    nOffC = nOffC+nB*nB
  end do
  call mma_deallocate(S1)
  call mma_deallocate(Occ)
  call mma_deallocate(CMO)

  if (Do_Cos) then
    Def = sum(Def_Irr(1:nSym))
    write(u6,'(/,6X,A,F9.4,A,I4,A,F6.3,A)') 'Deficit of the active space',Def,' of',nAct,' orbitals (threshold',ThrAct,')'
    Lost = (Def > ThrAct)
    if (Lost) then
      call WarningMessage(2,'The active space of the starting orbitals was lost.')
      write(u6,'(6X,A)') 'Active orbitals of the start have left the active space (see the weights above).'
      write(u6,'(6X,A)') 'Choose the active orbitals of the start anew (e.g. with a larger basis or a smaller step),'
      if (ActStop) then
        write(u6,'(6X,A)') 'or accept the result with EXPERT or a larger ACTThreshold. RASSCF stops with an error.'
      else
        write(u6,'(6X,A)') 'or raise ACTThreshold. EXPERT is given: RASSCF continues.'
      end if
    end if
  end if

  ! the occupation check: all active orbitals if the start has occupation
  ! numbers, else only the flagged ones
  nFlag = count(Flags(:) /= 0)
  if (Do_Match .or. (nFlag > 0)) then
    write(u6,*)
    write(u6,'(6X,A,F7.4,A)') 'Active natural orbitals: occupation numbers (flagged: above 2-t or below t, t =',ThrOcc,')'
    if (Do_Match) then
      write(u6,'(6X,A)') '(start: the occupation of the orbital in the starting active density; main starting orbital)'
      write(u6,'(6X,A)') 'Irrep  orbital  occupation  start occupation   main start orbital (overlap^2)'
    else
      write(u6,'(6X,A)') 'Irrep  orbital  occupation'
    end if
    do i=1,nAct
      if ((.not. Do_Match) .and. (Flags(i) == 0)) cycle
      if (Do_Match) then
        write(u6,'(6X,I5,I9,F12.4,F18.4,I13,A,F6.3,A,2X,A)') iSym_A(i),iOrb_A(i),Occ_Fin(i),Occ_Match(i),iMatch(i),' (',Ovl(i), &
                                                              ')',trim(Flag_Text(Flags(i),Contrib(i)))
      else
        write(u6,'(6X,I5,I9,F12.4,2X,A)') iSym_A(i),iOrb_A(i),Occ_Fin(i),trim(Flag_Text(Flags(i),Contrib(i)))
      end if
    end do
    if (nFlag > 0) then
      call WarningMessage(1,'Active natural orbitals with occupation numbers close to 0 or 2.')
      write(u6,'(6X,A)') 'Such an orbital adds hardly any correlation: it may have left the correlated space (e.g. a'
      write(u6,'(6X,A)') 'lone pair that became a pure atomic orbital). Check the active space; near a dissociation'
      write(u6,'(6X,A)') 'limit this can be intended. (Warning only; OCCThreshold 0 switches the check off.)'
    end if
  end if

#ifdef _HDF5_
  if (Have_H5) then
    if (Do_Cos) then
      call mh5_init_attr(H5Id,'ACTIVE_SPACE_DEFICIT',Def)
      call mh5_init_attr(H5Id,'ACTIVE_SPACE_DEFICIT_IRREP',1,[nSym],Def_Irr)
      call mh5_init_attr(H5Id,'ACTIVE_SPACE_COSINES',1,[nAct],Cos_All)
    end if
    if (ThrOcc > Zero) call mh5_init_attr(H5Id,'ACTIVE_SPACE_OCC_FLAGS',1,[nAct],Flags)
    if (Do_Match) call mh5_init_attr(H5Id,'ACTIVE_SPACE_OCC_START',1,[nAct],Occ_Match)
  end if
#endif
  call mma_deallocate(Contrib)
  call mma_deallocate(Ovl)
  call mma_deallocate(Occ_Match)
  call mma_deallocate(Occ_Fin)
  call mma_deallocate(iMatch)
  call mma_deallocate(iOrb_A)
  call mma_deallocate(iSym_A)
  call mma_deallocate(Flags)
  call mma_deallocate(Cos_All)
  call Free_Start()

end subroutine Act_Check_Final

!> The text after a flagged orbital: nearly full or nearly empty, and its
!> largest atomic contributions.
function Flag_Text(iFlag,Contrib) result(Txt)
  integer(kind=iwp), intent(in) :: iFlag
  character(len=*), intent(in) :: Contrib
  character(len=80) :: Txt
  select case (iFlag)
    case (1)
      Txt = '<-- nearly doubly occupied  '//Contrib
    case (-1)
      Txt = '<-- nearly empty  '//Contrib
    case default
      Txt = ''
  end select
end function Flag_Text

!> The largest atomic contributions (Mulliken gross populations, at most
!> three of at least 0.05) of an orbital C, SC = S C, Names: the basis
!> function labels (centre in the first LenIn characters).
subroutine Atom_Contrib(nB,C,SC,Names,Txt)

  use Molcas, only: LenIn
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: nB
  real(kind=wp), intent(in) :: C(nB), SC(nB)
  character(len=*), intent(in) :: Names(nB)
  character(len=*), intent(out) :: Txt
  integer(kind=iwp) :: i, iBest, iCen, l, nCen
  real(kind=wp), allocatable :: Pop(:)
  character(len=LenIn), allocatable :: Cen(:)
  character(len=16) :: Item

  call mma_allocate(Pop,nB,Label='Pop')
  call mma_allocate(Cen,nB,Label='Cen')
  nCen = 0
  do i=1,nB
    iCen = 0
    do l=1,nCen
      if (Cen(l) == Names(i)(1:LenIn)) then
        iCen = l
        exit
      end if
    end do
    if (iCen == 0) then
      nCen = nCen+1
      iCen = nCen
      Cen(iCen) = Names(i)(1:LenIn)
      Pop(iCen) = Zero
    end if
    Pop(iCen) = Pop(iCen)+C(i)*SC(i)
  end do
  Txt = ''
  do l=1,3
    if (nCen == 0) exit
    iBest = maxloc(Pop(1:nCen),1)
    if (Pop(iBest) < 0.05_wp) exit
    write(Item,'(A,1X,F4.2)') trim(Cen(iBest)),Pop(iBest)
    Txt = trim(Txt)//' '//trim(Item)
    Pop(iBest) = -huge(One)
  end do
  Txt = adjustl(Txt)
  call mma_deallocate(Cen)
  call mma_deallocate(Pop)

end subroutine Atom_Contrib

!> Free the copy of the starting orbitals (also when RASSCF ends before
!> the check, e.g. FCIDUMP only).
subroutine Act_Check_Free()
  call Free_Start()
end subroutine Act_Check_Free

subroutine Free_Start()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(CMO_Start,safe='*')
  call mma_deallocate(Occ_Start,safe='*')
  Have_Start = .false.
  Have_Occ_Start = .false.
  Do_Cos = .false.
end subroutine Free_Start

!> The singular values of the n x n matrix A, descending,
!> as square roots of the eigenvalues of A^T A. (dsyev instead of
!> dgesvd: dgesvd without singular vectors goes through dlasq2, whose
!> IEEE probe ieeeck traps in builds with floating-point exceptions.)
subroutine SVD_Values(n,A,SV)

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: A(n,n)
  real(kind=wp), intent(out) :: SV(n)
  integer(kind=iwp) :: Info, lWork
  real(kind=wp) :: Wrk1(1)
  real(kind=wp), allocatable :: ATA(:,:), Wrk(:)

  call mma_allocate(ATA,n,n,Label='ATA')
  call DGEMM_('T','N',n,n,n,One,A,n,A,n,Zero,ATA,n)
  call dsyev_('N','L',n,ATA,n,SV,Wrk1,-1,Info)
  lWork = max(1,int(Wrk1(1)))
  call mma_allocate(Wrk,lWork,Label='Wrk')
  call dsyev_('N','L',n,ATA,n,SV,Wrk,lWork,Info)
  call mma_deallocate(Wrk)
  call mma_deallocate(ATA)
  if (Info /= 0) then
    ! (no false alarm if the diagonalization fails)
    SV(:) = One
  else
    ! ascending eigenvalues -> descending singular values
    SV(:) = sqrt(max(SV(n:1:-1),Zero))
  end if

end subroutine SVD_Values

end module Active_Space_Check
