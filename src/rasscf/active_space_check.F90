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
module Active_Space_Check

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
private

! the threshold on the total deficit (ACTThreshold) and whether a lost
! active space stops the calculation (not with EXPERT)
real(kind=wp) :: ThrAct = 1.0_wp
logical(kind=iwp) :: ActStop = .true.
! copy of the starting orbitals (Act_Check_Start), the check applies
real(kind=wp), allocatable :: CMO_Start(:)
logical(kind=iwp) :: Have_Start = .false.
! the starting orbitals hold the active space of an earlier CASSCF/RASSCF
logical(kind=iwp) :: Restart = .false.
! active orbitals per irrep of the start file (type indices; -1: unknown)
integer(kind=iwp) :: nAct_File(8) = -1

public :: Act_Check_Final, Act_Check_Free, Act_Check_Start, Act_In_Loop, Act_Same_WF, Act_Set_Restart, ActStop, ThrAct

contains

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

!> Keep the orbitals the macro iterations start from (only in a restart).
subroutine Act_Check_Start(nCMO,CMO)

  use Orbital_Spaces, only: nAsh, nSym
  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nCMO
  real(kind=wp), intent(in) :: CMO(nCMO)

  Have_Start = .false.
  if ((ThrAct <= Zero) .or. (.not. Restart)) return
  ! a start file with another active space is not a restart of this one
  if (nAct_File(1) >= 0) then
    if (any(nAct_File(1:nSym) /= nAsh(1:nSym))) then
      write(u6,'(6X,A)') 'Active-space check skipped: the starting orbitals have another number of active orbitals'
      return
    end if
  end if
  call mma_allocate(CMO_Start,nCMO,Label='CMO_Start')
  CMO_Start(:) = CMO(:)
  Have_Start = .true.

end subroutine Act_Check_Start

!> Compare the converged active space (CMO) with the start; print the
!> weights outside of the start orbitals, the deficit and the cosines per
!> irrep, store the total deficit, the deficits per irrep and the cosines
!> on the HDF5 file (if open) and tell whether the active space was lost.
subroutine Act_Check_Final(nCMO,CMO,Have_H5,H5Id,Lost)

  use OneDat, only: sNoNuc, sNoOri
  use Index_Functions, only: nTri_Elem
  use Orbital_Spaces, only: nAsh, nBas, nFro, nIsh, nSym, nTot1
  use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _HDF5_
  use mh5, only: mh5_init_attr
#endif

  integer(kind=iwp), intent(in) :: nCMO, H5Id
  real(kind=wp), intent(in) :: CMO(nCMO)
  logical(kind=iwp), intent(in) :: Have_H5
  logical(kind=iwp), intent(out) :: Lost
  integer(kind=iwp) :: i, iC, iComp, iO, iOpt, iRc, iS, iSyLbl, j, k, nA, nAct, nB, nOff1, nOffC
  real(kind=wp) :: Def, Def_Irr(8), W
  real(kind=wp), allocatable :: CIn(:,:), COut(:,:), Cos_All(:), M(:,:), S1(:), SC(:,:), Sq(:,:), SV(:)
  character(len=8) :: Label
#include "macros.fh"

  Lost = .false.
#ifndef _HDF5_
  unused_var(H5Id)
  unused_var(Have_H5)
#endif
  if (.not. Have_Start) return
  nAct = sum(nAsh(1:nSym))
  if (nAct == 0) then
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
  if (iRc /= 0) then
    write(u6,'(/,6X,A)') 'Active-space check skipped: no AO overlap on ONEINT'
    call mma_deallocate(S1)
    call Free_Start()
    return
  end if

  call mma_allocate(Cos_All,nAct,Label='Cos_All')
  write(u6,*)
  write(u6,'(6X,A)') 'Active space: overlap of the converged active orbitals with the starting ones'
  write(u6,'(6X,A)') '(deficit = orbitals'' worth of the start that left the active space; 0 = unchanged)'
  nOff1 = 0
  nOffC = 0
  iC = 0
  Def_Irr(:) = Zero
  do iS=1,nSym
    nB = nBas(iS)
    nA = nAsh(iS)
    if (nA > 0) then
      call mma_allocate(Sq,nB,nB,Label='Sq')
      call Square(S1(nOff1+1),Sq,1,nB,nB)
      call mma_allocate(CIn,nB,nA,Label='CIn')
      call mma_allocate(COut,nB,nA,Label='COut')
      iO = nOffC+(nFro(iS)+nIsh(iS))*nB
      do j=1,nA
        CIn(:,j) = CMO_Start(iO+(j-1)*nB+1:iO+j*nB)
        COut(:,j) = CMO(iO+(j-1)*nB+1:iO+j*nB)
      end do
      call mma_allocate(SC,nB,nA,Label='SC')
      call mma_allocate(M,nA,nA,Label='M')
      call DGEMM_('N','N',nB,nA,nB,One,Sq,nB,COut,nB,Zero,SC,nB)
      call DGEMM_('T','N',nA,nA,nB,One,CIn,nB,SC,nB,Zero,M,nA)
      ! weights outside the converged active space (before the SVD
      ! overwrites M)
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
      call mma_allocate(SV,nA,Label='SV')
      call SVD_Values(nA,M,SV)
      write(u6,'(6X,A,F9.4)') 'deficit',Def_Irr(iS)
      write(u6,'(6X,A)') 'cosines of the principal angles:'
      write(u6,'(6X,8F9.4)') (SV(k),k=1,nA)
      Cos_All(iC+1:iC+nA) = SV(:)
      iC = iC+nA
      call mma_deallocate(SV)
      call mma_deallocate(M)
      call mma_deallocate(SC)
      call mma_deallocate(COut)
      call mma_deallocate(CIn)
      call mma_deallocate(Sq)
    end if
    nOff1 = nOff1+nTri_Elem(nB)
    nOffC = nOffC+nB*nB
  end do
  call mma_deallocate(S1)

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
#ifdef _HDF5_
  if (Have_H5) then
    call mh5_init_attr(H5Id,'ACTIVE_SPACE_DEFICIT',Def)
    call mh5_init_attr(H5Id,'ACTIVE_SPACE_DEFICIT_IRREP',1,[nSym],Def_Irr)
    call mh5_init_attr(H5Id,'ACTIVE_SPACE_COSINES',1,[nAct],Cos_All)
  end if
#endif
  call mma_deallocate(Cos_All)
  call Free_Start()

end subroutine Act_Check_Final

!> Free the copy of the starting orbitals (also when RASSCF ends before
!> the check, e.g. FCIDUMP only).
subroutine Act_Check_Free()
  call Free_Start()
end subroutine Act_Check_Free

subroutine Free_Start()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(CMO_Start,safe='*')
  Have_Start = .false.
end subroutine Free_Start

!> The singular values of the n x n matrix A (overwritten), descending.
subroutine SVD_Values(n,A,SV)

  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(inout) :: A(n,n)
  real(kind=wp), intent(out) :: SV(n)
  integer(kind=iwp) :: Info, lWork
  real(kind=wp) :: Dum(1,1), Wrk1(1)
  real(kind=wp), allocatable :: Wrk(:)

  call dgesvd_('N','N',n,n,A,n,SV,Dum,1,Dum,1,Wrk1,-1,Info)
  lWork = max(1,int(Wrk1(1)))
  call mma_allocate(Wrk,lWork,Label='Wrk')
  call dgesvd_('N','N',n,n,A,n,SV,Dum,1,Dum,1,Wrk,lWork,Info)
  call mma_deallocate(Wrk)
  ! (no false alarm if the SVD fails)
  if (Info /= 0) SV(:) = One

end subroutine SVD_Values

end module Active_Space_Check
