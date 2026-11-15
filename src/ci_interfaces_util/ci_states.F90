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

!> The gateway between the calling programs and the CI environments for the
!> CI states (split-graph walks, SGUGA): setup with explicit arguments,
!> operations by state index. The programs use only this module and CI_Info
!> (the copy of the state information); they never use the environments.
module CI_States

use sguga, only: CIS, EXS, SG_Checks, SG_Epq_Psi, SG_Free, SG_Init, SG_Init_Light, SGS, StepVec
use SG_Lucia_Bridge, only: SG_ReOrd, SG_ReOrd_Orig
use CI_Info, only: CI_St
use CI_Verify_Options, only: CI_Verify
use CSF_Graph, only: SG_Walk_Next
use Molcas, only: MxLev, MxSym
use Definitions, only: wp, iwp

implicit none
private

public :: CI_CITRA, CI_Epq_Psi, CI_Free_State, CI_Occ_Diag, CI_Print_Graph, CI_PrWF, CI_ReOrd, CI_ReOrd_Orig, CI_Setup_SG, &
          CI_StepVec, CI_SYG2SGU, CI_Walk_Next

contains

!> Set up the split-graph CI space of state iState (SGUGA): nSym, nActEl,
!> iSpin; the RAS partitions nRas(iSym,iPrt) with the electron limits
!> nRasEl(iPrt); Level/L2Act (orbital <-> level), nLev levels with
!> symmetries NSM. Mode 'FULL' (graph, walks and coupling tables), 'GRAPH'
!> (graph and walks), 'NONE' (no graph; with LSM_Single the CI space is
!> the single function of symmetry LSM_Single). Fills CI_St(iState).
subroutine CI_Setup_SG(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,Level,L2Act,nLev,NSM,Mode,LSM_Single)
  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin, nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt), Level(MxLev), &
                                   L2Act(MxLev), nLev, NSM(MxLev)
  character(len=*), intent(in) :: Mode
  integer(kind=iwp), intent(in), optional :: LSM_Single

  SG_Checks = CI_Verify('ANY')
  select case (Mode)
    case ('FULL')
      call SG_Init(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM)
    case ('GRAPH')
      call SG_Init_Light(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM)
    case ('NONE')
      call SG_Init_Light(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM, &
                         Do_MkSGUGA=.false.)
      if (present(LSM_Single)) CIS(iState)%NCSF(LSM_Single) = 1
    case default
      call Abend()
  end select
  call Fill_Info(iState)
end subroutine CI_Setup_SG

!> The copy of the state information for the programs.
subroutine Fill_Info(iState)
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: n
  CI_St(iState)%Active = .true.
  CI_St(iState)%nSym = SGS(iState)%nSym
  CI_St(iState)%nActEl = SGS(iState)%nActEl
  CI_St(iState)%iSpin = SGS(iState)%iSpin
  CI_St(iState)%nLev = SGS(iState)%nLev
  CI_St(iState)%MidLev = SGS(iState)%MidLev
  CI_St(iState)%IFRAS = SGS(iState)%IFRAS
  CI_St(iState)%nRsPrt = SGS(iState)%nRsPrt
  CI_St(iState)%nCSF(:) = 0
  if (allocated(CIS(iState)%NCSF)) then
    n = min(size(CIS(iState)%NCSF),MxSym)
    CI_St(iState)%nCSF(1:n) = CIS(iState)%NCSF(1:n)
  end if
  CI_St(iState)%Level(:) = SGS(iState)%Level(:)
  CI_St(iState)%L2Act(:) = SGS(iState)%L2Act(:)
  CI_St(iState)%ISm(:) = 0
  if (allocated(SGS(iState)%ISm)) CI_St(iState)%ISm(1:SGS(iState)%nLev) = SGS(iState)%ISm(1:SGS(iState)%nLev)
end subroutine Fill_Info

!> Free the CI environment of state iState and its information.
subroutine CI_Free_State(iState)
  integer(kind=iwp), intent(in) :: iState
  call SG_Free(iState)
  CI_St(iState)%Active = .false.
  CI_St(iState)%nCSF(:) = 0
end subroutine CI_Free_State

!> SGM = SGM + CPQ E_pq CI for state iState (CI of symmetry ISYCI).
subroutine CI_Epq_Psi(iState,IP,IQ,CPQ,ISYCI,CI,SGM)
  integer(kind=iwp), intent(in) :: iState, IP, IQ, ISYCI
  real(kind=wp), intent(in) :: CPQ, CI(*)
  real(kind=wp), intent(inout) :: SGM(*)
  call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState),IP,IQ,CPQ,ISYCI,CI,SGM)
end subroutine CI_Epq_Psi

!> The CSFs of symmetry LSM of state iState one by one in the order of the
!> CI vectors: it(1:4) = (MV, ISYUP, IDWN, IUP) of the next CSF (start with
!> it = 1; it(1) = 0 after the last one); ICS(1:nLev) its step vector
!> (0 empty, 1 up, 2 down, 3 doubly occupied). Before the call, it(:) is
!> the split-graph position (midvertex, upper symmetry, lower and upper
!> walk) of the CSF returned.
subroutine CI_Walk_Next(iState,LSM,it,ICS)
  integer(kind=iwp), intent(in) :: iState, LSM
  integer(kind=iwp), intent(inout) :: it(4)
  integer(kind=iwp), intent(out) :: ICS(*)
  call SG_Walk_Next(SGS(iState),CIS(iState),LSM,it(1),it(2),it(3),it(4),ICS)
end subroutine CI_Walk_Next

!> DIAG(I) = sum_lev n_lev(I) ETA(lev) for the CSFs I of symmetry ISYCI
!> (occupation numbers n of the levels).
subroutine CI_Occ_Diag(iState,ISYCI,ETA,DIAG,nDiag)
  integer(kind=iwp), intent(in) :: iState, ISYCI, nDiag
  real(kind=wp), intent(in) :: ETA(*)
  real(kind=wp), intent(out) :: DIAG(nDiag)
  call SG_Occ_Diag(iState,ISYCI,ETA,DIAG,nDiag)
end subroutine CI_Occ_Diag

!> Print the CI vector CI of symmetry ISYCI of state iState (SG_PRWF).
subroutine CI_PrWF(iState,ISYCI,CITHR,iSpin,CI,lCI,KeyPRSD,LuVecDet)
  integer(kind=iwp), intent(in) :: iState, ISYCI, iSpin, lCI, LuVecDet
  logical(kind=iwp), intent(in) :: KeyPRSD
  real(kind=wp), intent(in) :: CITHR, CI(lCI)
  call SG_PRWF(iState,ISYCI,CITHR,iSpin,CI,lCI,KeyPRSD,LuVecDet)
end subroutine CI_PrWF

!> Print the split graph of state iState.
subroutine CI_Print_Graph(iState)
  integer(kind=iwp), intent(in) :: iState
  call SG_Print(iState)
end subroutine CI_Print_Graph

!> The CI transformation for an orbital transformation (SG_CITRA).
subroutine CI_CITRA(WFTP,iState,LSM,NTRA,TRA,NCO,CI,NOSH,NISH,NASH)
  character(len=8), intent(in) :: WFTP
  integer(kind=iwp), intent(in) :: iState, LSM, NCO, NTRA, NOSH(MxSym), NISH(MxSym), NASH(MxSym)
  real(kind=wp), intent(in) :: TRA(NTRA)
  real(kind=wp), intent(inout) :: CI(NCO)
  call SG_CITRA(WFTP,iState,LSM,NTRA,TRA,NCO,CI,NOSH,NISH,NASH)
end subroutine CI_CITRA

!> Reorder a CI vector between the Lucia and the split-graph order.
subroutine CI_ReOrd(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
  integer(kind=iwp), intent(in) :: iState, IREFSM, IMODE, nConf
  real(kind=wp), intent(in) :: CIOLD(nConf)
  real(kind=wp), intent(out) :: CINEW(nConf)
  call SG_ReOrd(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
end subroutine CI_ReOrd

!> As CI_ReOrd, from Lucia's configuration tables (MCLR).
subroutine CI_ReOrd_Orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
  integer(kind=iwp), intent(in) :: iState, IREFSM, IMODE, nConf
  real(kind=wp), intent(in) :: CIOLD(nConf)
  real(kind=wp), intent(out) :: CINEW(nConf)
  call SG_ReOrd_Orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
end subroutine CI_ReOrd_Orig

!> The walk (step vector) of a configuration with the closed shells ICL,
!> open shells IOP and their spin coupling ISPIN (STEPVEC).
subroutine CI_StepVec(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)
  integer(kind=iwp), intent(in) :: NCL, NOP, NORB, ICL(NCL), IOP(NOP), ISPIN(NOP)
  integer(kind=iwp), intent(out) :: IWALK(NORB)
  call STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)
end subroutine CI_StepVec

!> Transform a CI vector between the symmetric-group (RASSI) and the
!> split-graph basis of state iState (SYG2SGU).
subroutine CI_SYG2SGU(IMODE,iState,LSYM,ICNFTAB,ISPNTAB,CIOLD,CINEW,iPrint)
  integer(kind=iwp), intent(in) :: IMODE, iState, LSYM, ICNFTAB(*), ISPNTAB(*), iPrint
  real(kind=wp), intent(in) :: CIOLD(*)
  real(kind=wp), intent(out) :: CINEW(*)
  call SYG2SGU(IMODE,iState,LSYM,ICNFTAB,ISPNTAB,CIOLD,CINEW,iPrint)
end subroutine CI_SYG2SGU

end module CI_States
