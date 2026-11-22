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
#include "macros.fh"

module CI_States

use sguga, only: CIS, EXS, SG_Checks, SG_Epq_Psi, SG_Free, SG_Init, SG_Init_Light, SG_New_Tag, SG_No_Reuse, SGS, StepVec
use SG_Lucia_Bridge, only: SG_ReOrd, SG_ReOrd_Orig
use SGUGA_Sigma, only: SG_Occ_Diag
use SGUGA_Print, only: SG_Print, SG_PRWF
use SGUGA_Transform, only: SG_CITRA, SYG2SGU
use CI_Info, only: CI_St, Space_CI
use CI_Space_Def, only: Space_Free, Space_Pack, Space_Rec_Len, Space_Unpack
use CI_Verify_Options, only: CI_Verify
use CSF_Graph, only: SG_Walk_Next
use Molcas, only: MxLev, MxSym
use Definitions, only: wp, iwp

implicit none
private

public :: CI_CITRA, CI_Epq_Psi, CI_Free_State, CI_Occ_Diag, CI_Print_Graph, CI_PrWF, CI_ReOrd, CI_ReOrd_Orig, CI_Setup_RAS, &
          CI_Setup_SG, CI_StepVec, CI_SYG2SGU, CI_Walk_Next, RAS_Partitions
public :: CI_New_Tag, CI_No_GAS, CI_No_Reuse, CI_Space_Free, CI_Space_From_Rec, CI_Space_Rd_H5, CI_Space_Rd_JobIph, &
          CI_Space_Rec_Len, CI_Space_To_Rec

! length of the CI-space record of JobIph (IADR15(19)) and HDF5 (CI_SPACE)
integer(kind=iwp), parameter :: CI_Space_Rec_Len = Space_Rec_Len
!> Tag of a CI vector for CI_Epq_Psi: a caller that applies several E_pq to
!> the same vector passes the same tag (a new one, from CI_New_Tag, whenever
!> the vector changes), so that intermediates can be reused without a pass
!> over the vector. CI_No_Reuse: the vector is used only once.
integer(kind=iwp), parameter :: CI_No_Reuse = SG_No_Reuse

contains

!> The partitions of a RAS space (RAS1, RAS2, RAS3 with nRs1, nRs2, nRs3
!> orbitals per symmetry, at most nHole1 holes in RAS1 and nElec3 electrons
!> in RAS3) in the form of CI_Setup_SG: nPrt partitions with nRas(iSym,k)
!> orbitals and at least nRasEl(k) electrons in the partitions 1..k. With
!> holes or electrons restricted, or All_Parts, the three partitions;
!> otherwise (CAS) RAS2 as the only partition.
subroutine RAS_Partitions(nSym,nActEl,nRs1,nRs2,nRs3,nHole1,nElec3,All_Parts,nPrt,nRas,nRasEl)
  integer(kind=iwp), intent(in) :: nSym, nActEl, nRs1(nSym), nRs2(nSym), nRs3(nSym), nHole1, nElec3
  logical(kind=iwp), intent(in) :: All_Parts
  integer(kind=iwp), intent(out) :: nPrt, nRas(MxSym,3), nRasEl(3)
  nRas(:,:) = 0
  nRasEl(:) = 0
  if (All_Parts .or. (nHole1+nElec3 /= 0)) then
    nPrt = 3
    nRas(1:nSym,1) = nRs1(1:nSym)
    nRas(1:nSym,2) = nRs2(1:nSym)
    nRas(1:nSym,3) = nRs3(1:nSym)
    nRasEl(1) = 2*sum(nRs1(1:nSym))-nHole1
    nRasEl(2) = nActEl-nElec3
    nRasEl(3) = nActEl
  else
    nPrt = 1
    nRas(1:nSym,1) = nRs2(1:nSym)
    nRasEl(1) = nActEl
  end if
end subroutine RAS_Partitions

!> CI_Setup_SG for a RAS space (see RAS_Partitions): Level(1:nLev) and
!> L2Act(1:nLev) the orbital <-> level maps, NSM(1:nLev) the symmetries of
!> the levels; Mode and LSM_Single as in CI_Setup_SG.
subroutine CI_Setup_RAS(iState,nSym,nActEl,iSpin,LSym,nRs1,nRs2,nRs3,nHole1,nElec3,All_Parts,nLev,Level,L2Act,NSM,Mode,LSM_Single)
  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin, LSym, nRs1(nSym), nRs2(nSym), nRs3(nSym), nHole1, nElec3, nLev, &
                                   Level(nLev), L2Act(nLev), NSM(nLev)
  logical(kind=iwp), intent(in) :: All_Parts
  character(len=*), intent(in) :: Mode
  integer(kind=iwp), intent(in), optional :: LSM_Single
  integer(kind=iwp) :: nPrt, nRas(MxSym,3), nRasEl(3), LevelX(MxLev), L2ActX(MxLev), NSMX(MxLev), iq
  call RAS_Partitions(nSym,nActEl,nRs1,nRs2,nRs3,nHole1,nElec3,All_Parts,nPrt,nRas,nRasEl)
  LevelX(:) = [(iq,iq=1,MxLev)]
  L2ActX(:) = [(iq,iq=1,MxLev)]
  NSMX(:) = 0
  LevelX(1:nLev) = Level(1:nLev)
  L2ActX(1:nLev) = L2Act(1:nLev)
  NSMX(1:nLev) = NSM(1:nLev)
  if (present(LSM_Single)) then
    call CI_Setup_SG(iState,nSym,nActEl,iSpin,LSym,nRas(:,1:nPrt),nRasEl(1:nPrt),nPrt,LevelX,L2ActX,nLev,NSMX,Mode, &
                     LSM_Single=LSM_Single)
  else
    call CI_Setup_SG(iState,nSym,nActEl,iSpin,LSym,nRas(:,1:nPrt),nRasEl(1:nPrt),nPrt,LevelX,L2ActX,nLev,NSMX,Mode)
  end if
end subroutine CI_Setup_RAS

!> Set up the split-graph CI space of state iState (SGUGA): nSym, nActEl,
!> iSpin, LSym the symmetry of the CI state; the RAS partitions nRas(iSym,iPrt) with the electron limits
!> nRasEl(iPrt); Level/L2Act (orbital <-> level), nLev levels with
!> symmetries NSM. Mode 'FULL' (graph, walks and coupling tables), 'GRAPH'
!> (graph and walks), 'NONE' (no graph; with LSM_Single the CI space is
!> the single function of symmetry LSM_Single). Fills CI_St(iState).
subroutine CI_Setup_SG(iState,nSym,nActEl,iSpin,LSym,nRas,nRasEl,nRsPrt,Level,L2Act,nLev,NSM,Mode,LSM_Single,nRasMx)
  integer(kind=iwp), intent(in) :: iState, nSym, nActEl, iSpin, LSym, nRsPrt, nRas(MxSym,nRsPrt), nRasEl(nRsPrt), Level(MxLev), &
                                   L2Act(MxLev), nLev, NSM(MxLev)
  character(len=*), intent(in) :: Mode
  ! LSM_Single: see above; nRasMx: GAS, upper bounds on the electrons in the
  ! partitions 1..k (nRasEl the lower ones)
  integer(kind=iwp), intent(in), optional :: LSM_Single, nRasMx(nRsPrt)

  SG_Checks = CI_Verify('ANY')
  select case (Mode)
    case ('FULL')
      call SG_Init(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM,xnRasMx=nRasMx)
    case ('GRAPH')
      call SG_Init_Light(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM, &
                         xnRasMx=nRasMx)
    case ('NONE')
      call SG_Init_Light(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=NSM, &
                         Do_MkSGUGA=.false.,xnRasMx=nRasMx)
      if (present(LSM_Single)) CIS(iState)%NCSF(LSM_Single) = 1
    case default
      call Abend()
  end select
  call Fill_Info(iState)
  ! the CI state itself, also without a graph (Mode 'NONE' marks the graph
  ! with SGS%nActEl = SGS%iSpin = 0)
  CI_St(iState)%nActEl = nActEl
  CI_St(iState)%iSpin = iSpin
  CI_St(iState)%STSym = LSym
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
  CI_St(iState)%STSym = 0
end subroutine CI_Free_State

!> SGM = SGM + CPQ E_pq CI for state iState (CI of symmetry ISYCI); CI_Tag:
!> see CI_New_Tag (without it the vector is recognized by its contents).
subroutine CI_Epq_Psi(iState,IP,IQ,CPQ,ISYCI,CI,SGM,CI_Tag)
  integer(kind=iwp), intent(in) :: iState, IP, IQ, ISYCI
  real(kind=wp), intent(in) :: CPQ, CI(*)
  real(kind=wp), intent(inout) :: SGM(*)
  integer(kind=iwp), intent(in), optional :: CI_Tag
  call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState),IP,IQ,CPQ,ISYCI,CI,SGM,CI_Tag)
end subroutine CI_Epq_Psi

!> A new tag for a CI vector (see CI_Epq_Psi).
function CI_New_Tag() result(iTag)
  integer(kind=iwp) :: iTag
  iTag = SG_New_Tag()
end function CI_New_Tag

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

!> The CI space of the wave function (CI_Info's Space_CI) as the record of
!> JobIph/HDF5; GAS tells whether it came from a GAS input. Have is false
!> if there is no space (e.g. DMRG, NECI: no CI of the CI interface).
subroutine CI_Space_To_Rec(GAS,n,Rec,Have)
  logical(kind=iwp), intent(in) :: GAS
  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(out) :: Rec(n)
  logical(kind=iwp), intent(out) :: Have
  Have = (Space_CI%nPart > 0)
  Rec(:) = 0
  if (Have) call Space_Pack(Space_CI,GAS,n,Rec)
end subroutine CI_Space_To_Rec

!> CI_Info's Space_CI from the record of JobIph/HDF5; GAS tells whether the
!> wave function came from a GAS input.
subroutine CI_Space_From_Rec(n,Rec,GAS)
  integer(kind=iwp), intent(in) :: n, Rec(n)
  logical(kind=iwp), intent(out) :: GAS
  call Space_Unpack(n,Rec,Space_CI,GAS)
end subroutine CI_Space_From_Rec

!> The CI-space record of a JobIph at the disk address iAdr (IADR15(19));
!> iAdr <= 0 (files without the record): Found = .false. and no space, the
!> reader keeps the RAS description of the header.
subroutine CI_Space_Rd_JobIph(Lu,iAdr,Found,GAS)
  integer(kind=iwp), intent(in) :: Lu, iAdr
  logical(kind=iwp), intent(out) :: Found, GAS
  integer(kind=iwp) :: iDisk, Rec(CI_Space_Rec_Len)
  GAS = .false.
  Found = (iAdr > 0)
  if (.not. Found) then
    call Space_Free(Space_CI)
    return
  end if
  iDisk = iAdr
  call iDaFile(Lu,2,Rec,CI_Space_Rec_Len,iDisk)
  call CI_Space_From_Rec(CI_Space_Rec_Len,Rec,GAS)
end subroutine CI_Space_Rd_JobIph

!> The CI-space record of an HDF5 wave-function file (attribute CI_SPACE);
!> files without it (older files, DMRG, NECI) give Found = .false.; their
!> GAS flag is then taken from CI_TYPE.
subroutine CI_Space_Rd_H5(id,Found,GAS)
# ifdef _HDF5_
  use mh5, only: mh5_exists_attr, mh5_fetch_attr
# endif
  integer(kind=iwp), intent(in) :: id
  logical(kind=iwp), intent(out) :: Found, GAS
# ifdef _HDF5_
  integer(kind=iwp) :: Rec(CI_Space_Rec_Len)
  character(len=3) :: CIType
# endif
  Found = .false.
  GAS = .false.
# ifdef _HDF5_
  if (mh5_exists_attr(id,'CI_SPACE')) then
    call mh5_fetch_attr(id,'CI_SPACE',Rec)
    call CI_Space_From_Rec(CI_Space_Rec_Len,Rec,GAS)
    Found = .true.
  else
    call Space_Free(Space_CI)
    if (mh5_exists_attr(id,'CI_TYPE')) then
      call mh5_fetch_attr(id,'CI_TYPE',CIType)
      GAS = (CIType == 'GAS')
    end if
  end if
# else
  unused_var(id)
# endif
end subroutine CI_Space_Rd_H5

!> Free CI_Info's Space_CI (read by CI_Space_Rd_*).
subroutine CI_Space_Free()
  call Space_Free(Space_CI)
end subroutine CI_Space_Free

!> Stop a program that cannot yet treat a GAS reference wave function.
subroutine CI_No_GAS(ProgName)
  use Definitions, only: u6
  character(len=*), intent(in) :: ProgName
  character(len=len(ProgName)) :: Name
  Name = ProgName
  call UpCase(Name)
  call WarningMessage(2,trim(Name)//': GAS reference wave functions are not supported yet.')
  write(u6,'(1X,A)') 'The reference wave function was optimized with GASSCF (CI space with general'
  write(u6,'(1X,A)') 'active-space partitions), which '//trim(Name)//' cannot yet use.'
  call Quit_OnUserError()
end subroutine CI_No_GAS

end module CI_States
