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

module CI_Interfaces
use sguga, only: SG_Free
use SG_Lucia_Bridge, only: SG_Lucia_Walk, SG_ReOrd, SG_Verify_Lucia_Walk
use Lucia_Interface, only: Lucia_Util
use faroald, only: faroald_full, my_norb, sigma_update, htu, gtuvx, ndeta, ndetb ,transition_one_pdm, one_pdm, two_pdm, &
                   fold_two_pdm, traci_faroald
use citrans, only: citrans_conf_dets, citrans_csf2ras, citrans_csf2sd, citrans_ncsf, citrans_ras2csf, citrans_sd2csf, citrans_sort
use faroald_ras, only: faroald_ras_free, faroald_ras_init, n_ras, ras_active, ras_sigma, ras_traci, ras_trans_one_pdm, &
                       ras_trans_two_pdm
use rasscf_global, only: DoFaro, DoSGEx, DoSGCI
use general_data, only: nAsh
use rasscf_global, only: NAC, nDet, nConf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use definitions, only: wp, iwp, u6
use sguga, only: CIS, SGS, EXS, SGCnfList, SG_Explicit_Block, SG_Free_Cnf_List
use SGUGA_Sigma, only: sg_h_psi, sg_h_psi_ras, sg_hdiag_fast, sg_hdiag_ref
use SGUGA_Densities, only: sg_one_pdm, sg_spin_pdm, sg_trans_pdm, sg_two_pdm
use general_data, only: iDoGAS
use CI_Verify_Options, only: CI_Verify
use CI_Space_Def, only: Space_Free, Space_From_GAS, Space_Max_Occ, Space_Print, Space_Remove, Space_Removals
use CI_Info, only: CI_St, Space_CI
use Molcas, only: MxGAS

#include "intent.fh"

implicit none
Private

Public :: Mk_Tra_Psi, Mk_H_Psi, Mk_H_Psi_Det, Mk_PDMS, Mk_T1DM, Mk_IHI, Mk_IHJ, Mk_Cnf_Layout
Public :: CI_Basis, CI_Basis_To_File, CI_Close, CI_Faro_nDet_RAS, CI_File_Index_To_Basis, CI_File_To_Basis, CI_Initialize, &
          CI_Lucia_Active, CI_Timer, CI_Verify, Faro_nTmp, SG_CI

! CI_Verify (module CI_Verify_Options) is re-exported, see there for the
! run-time verification options (MOLCAS_CI_VERIFY).

! Mk_IHJ is verified against Mk_H_Psi only once per run.
logical(kind=iwp) :: IHJ_Verified = .false.
! sg_h_psi_ras is compared with sg_h_psi (CAS) only once per run; for RAS
! the difference to the CI back end is printed once (checked every time).
logical(kind=iwp) :: HPsi_RAS_Verified = .false.
! The SGUGA CI transformation (Mk_Tra_Psi, SIGMA) is checked on every call,
! its maximum error printed once.
logical(kind=iwp) :: TraPsi_Verified = .false.
! SGEX given but SGUGA not available: warned once.
logical(kind=iwp) :: SGEx_Warned = .false.
! SGCI: the SGUGA sigma vector compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: HPsi_SG_Verified = .false.
! SGCI: the SGUGA configuration list of the CI symmetry (SG_Cnf_Kept)
type(SGCnfList) :: Cnf_SG
! SGCI: the SGUGA densities compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: PDM_SG_Verified = .false.
! SGCI: the SGUGA diagonal checked against sg_hdiag_ref (SIGMA), printed once.
logical(kind=iwp) :: IHI_SG_Verified = .false.
! SGCI: the SGUGA configuration-average diagonal checked against Lucia, printed once.
logical(kind=iwp) :: IHI_Avg_Verified = .false.
! SGCI: the SGUGA transition densities compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: T1DM_SG_Verified = .false.
! FARO: the Faroald steps verified against Lucia and SGUGA (SIGMA), printed once
logical(kind=iwp) :: HPsi_Faro_Verified = .false., PDM_Faro_Verified = .false., T1DM_Faro_Verified = .false., &
                     TraPsi_Faro_Verified = .false., IHI_Faro_SG_Verified = .false.
! Lucia has been initialized (not for SGCI, except as the SIGMA reference)
logical(kind=iwp) :: Lucia_Active = .true.
! FARO: Faroald works with the active orbitals in the order of the levels
! of the split graph (by RAS partition, within each by symmetry), so that
! its CSFs couple the spins in the order of the file basis. Faro_L2O(lev)
! is the active orbital (symmetry-blocked order of TU/TUVX) of level lev;
! the identity for CAS. The densities are in level order (as Lucia's and
! SGUGA's).
integer(kind=iwp), allocatable :: Faro_L2O(:)
! FARO with RAS: the string-driven sigma (faroald_ras) checked once against
! the sigma of the full determinant space (SIGMA)
logical(kind=iwp) :: HPsi_FaroRAS_Verified = .false.
! SIGMA: the spin density of the ASTA triplet operator checked, printed once
logical(kind=iwp) :: SD_ASTA_Verified = .false.
! HEXS/DEXS (RASSCF keywords, CI_Info Space_Inp): the configurations
! with the maximum occupation (HEXS) or at least the maximum - 1 (DEXS) in
! the listed RAS partitions are removed as in Lucia: their CSFs stay in the
! vectors, the CI works with P H P (zero sigma, zero diagonal there). The
! maximum occupation of partition g is maxa(g) + maxb(g) (Space_Max_Occ), the
! largest numbers of alpha and beta electrons in g in the determinants of
! the space (Lucia's rule for the determinant blocks, I_AM_SO_EXCITED).
! HEXS_Keep: the kept CSFs in the CI basis (SGCI, FARO; Lucia does it
! itself); HEXS_Keep_SG: in the split-graph basis (the SGUGA references).
logical(kind=iwp), allocatable :: HEXS_Keep(:), HEXS_Keep_SG(:)

! GAS: SGUGA evaluates products of excitation operators (H, 2-RDM, spin and
! transition densities, diagonal, explicit block) in the relaxed space of
! the state (CI_St(iState)%iRel, see SG_Rel_State), which holds all their
! intermediate states; Rel_Map(LSM)%M(I) is the index there of the CSF I of
! symmetry LSM (SG_Rel_Map).
type Rel_Map_Type
  integer(kind=iwp), allocatable :: M(:)
end type Rel_Map_Type
type(Rel_Map_Type) :: Rel_Map(8)

contains

!> SGCI: the whole CI (diagonal, explicit subspace, sigma vectors, density
!> matrices) is done by SGUGA in its own split-graph CSF basis, without
!> Lucia (RASSCF keyword SGCI, only where SGUGA covers the CI space).
function SG_CI() result(SGMode)
logical(kind=iwp) :: SGMode
SGMode = DoSGCI .and. SG_Avail()
end function SG_CI

!> FARO: the number of determinants of the RAS/GAS space of Faroald (all
!> symmetries; nela/nelb electrons, o(k) orbitals in partition k of np,
!> mn(k) to mx(k) electrons in the partitions 1..k), before any setup, for
!> the memory estimate of RASSCF.
function CI_Faro_nDet_RAS(nela,nelb,np,o,mn,mx) result(n)
use faroald_ras, only: ras_count_dets
integer(kind=iwp), intent(in) :: nela, nelb, np, o(np), mn(np), mx(np)
real(kind=wp) :: n
n = ras_count_dets(nela,nelb,np,o,mn,mx)
end function CI_Faro_nDet_RAS

!> FARO: the size of the determinant scratch of Faroald_Sigma (ctemp,
!> sigtemp): the full determinant matrix, not needed with RAS (faroald_ras).
function Faro_nTmp() result(n)
integer(kind=iwp) :: n
if (ras_active) then
  n = 1
else
  n = nDetA*nDetB
end if
end function Faro_nTmp

!> Lucia is initialized (CI_Initialize): its tables (determinant counts,
!> configuration types) may be used.
function CI_Lucia_Active() result(Active)
logical(kind=iwp) :: Active
Active = Lucia_Active
end function CI_Lucia_Active

!> The CSF basis of the CI vectors of RASSCF: 1 = Lucia, 2 = the split-graph
!> basis of SGUGA (SGCI; the order of the files), 3 = the configuration
!> order of Faroald (citrans; keyword FARO).
function CI_Basis() result(iBasis)
integer(kind=iwp) :: iBasis
if (SG_CI()) then
  iBasis = 2
else if (DoFaro) then
  iBasis = 3
else
  iBasis = 1
end if
end function CI_Basis

!> A CI vector from the basis of the CI back end to the split-graph order of
!> the JOBIPH/HDF5 files (and of SGUGA).
subroutine CI_Basis_To_File(iState,LSM,n,Cin,Cout)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: Cin(n)
real(kind=wp), intent(out) :: Cout(n)
select case (CI_Basis())
  case (1)
    call SG_Reord(iState,LSM,0,n,Cin,Cout)
  case (2)
    Cout(:) = Cin(:)
  case (3)
    call CITRANS_SORT('O',Cin,Cout)
end select
end subroutine CI_Basis_To_File

!> A CI vector from the split-graph order of the files to the basis of the
!> CI back end.
subroutine CI_File_To_Basis(iState,LSM,n,Cin,Cout)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: Cin(n)
real(kind=wp), intent(out) :: Cout(n)
select case (CI_Basis())
  case (1)
    call SG_Reord(iState,LSM,1,n,Cin,Cout)
  case (2)
    Cout(:) = Cin(:)
  case (3)
    call CITRANS_SORT('C',Cin,Cout)
end select
end subroutine CI_File_To_Basis

!> A CI vector from the basis of the CI back end to the Lucia basis (the
!> reference of MOLCAS_CI_VERIFY=SIGMA).
subroutine CI_Basis_To_Lucia(iState,LSM,n,Cin,Cout)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: Cin(n)
real(kind=wp), intent(out) :: Cout(n)
real(kind=wp), allocatable :: Tmp(:)
if (CI_Basis() == 1) then
  Cout(:) = Cin(:)
else
  call mma_allocate(Tmp,n,Label='Tmp_File')
  call CI_Basis_To_File(iState,LSM,n,Cin,Tmp)
  call SG_Reord(iState,LSM,1,n,Tmp,Cout)
  call mma_deallocate(Tmp)
end if
end subroutine CI_Basis_To_Lucia

!> A CI vector from the Lucia basis to the basis of the CI back end.
subroutine CI_Lucia_To_Basis(iState,LSM,n,Cin,Cout)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: Cin(n)
real(kind=wp), intent(out) :: Cout(n)
real(kind=wp), allocatable :: Tmp(:)
if (CI_Basis() == 1) then
  Cout(:) = Cin(:)
else
  call mma_allocate(Tmp,n,Label='Tmp_File')
  call SG_Reord(iState,LSM,0,n,Cin,Tmp)
  call CI_File_To_Basis(iState,LSM,n,Tmp,Cout)
  call mma_deallocate(Tmp)
end if
end subroutine CI_Lucia_To_Basis

!> Reference sigma vectors (CI development, MOLCAS_CI_VERIFY=SIGMA): Sout =
!> H Cin, both in the basis of the CI back end, computed by Lucia (Lab =
!> 'Lucia', the vector reordered to the Lucia basis and back) or by SGUGA
!> (Lab = 'SGUGA', through the split-graph order of the files).
subroutine Ref_Sigma(Lab,iState,LSM,n,Cin,Sout,nTU,TU,nTUVX,TUVX)
use lucia_data, only: IREFSM, Sigma_on_disk
use cands, only: ICSM, ISSM
character(len=*), intent(in) :: Lab
integer(kind=iwp), intent(in) :: iState, LSM, n, nTU, nTUVX
real(kind=wp), intent(in) :: Cin(n), TU(nTU), TUVX(nTUVX)
real(kind=wp), intent(out) :: Sout(n)
integer(kind=iwp) :: IREFSM_Keep, nD, nT
logical(kind=iwp) :: Sigma_Keep
real(kind=wp), allocatable :: C1(:), C2(:)

select case (Lab)
  case ('Lucia')
    nD = CI_St(iState)%nDet(LSM)
    nT = max(n,nD)
    call mma_allocate(C1,nT,Label='C1')
    call mma_allocate(C2,nT,Label='C2')
    IREFSM_Keep = IREFSM
    Sigma_Keep = Sigma_on_disk
    IREFSM = LSM
    ICSM = LSM
    ISSM = LSM
    C1(:) = Zero
    call CI_Basis_To_Lucia(iState,LSM,n,Cin,C1(1:n))
    C2(:) = Zero
    call csdtvc(C1,C2,1,LSM,0)
    C1(1:nD) = C2(1:nD)
    call Lucia_Util('Sigma',CI_Vector=C1(:),Sigma_Vector=C2(:),nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
    call CSDTVC(C1,C2,2,LSM,0)
    call CI_Lucia_To_Basis(iState,LSM,n,C1(1:n),Sout)
    IREFSM = IREFSM_Keep
    ICSM = IREFSM
    ISSM = IREFSM
    Sigma_on_disk = Sigma_Keep
  case ('SGUGA')
    call mma_allocate(C1,n,Label='C1')
    call mma_allocate(C2,n,Label='C2')
    call CI_Basis_To_File(iState,LSM,n,Cin,C1)
    call HEXS_Proj(HEXS_Keep_SG,C1)
    call SG_HPsi(iState,C1,n,LSM,C2,TUVX,nTUVX,TU,nTU,CI_St(iState)%nRsPrt /= 1)
    call HEXS_Proj(HEXS_Keep_SG,C2)
    call CI_File_To_Basis(iState,LSM,n,C2,Sout)
  case default
    write(u6,*) 'Ref_Sigma: unknown back end ',Lab
    call Abend()
end select
call mma_deallocate(C2)
call mma_deallocate(C1)
end subroutine Ref_Sigma

!> Report (once, Done) and check a maximal deviation of a verification.
subroutine Report_Dev(Text,Dev,Thr,Done)
character(len=*), intent(in) :: Text
real(kind=wp), intent(in) :: Dev, Thr
logical(kind=iwp), intent(inout) :: Done
if (.not. Done) write(u6,'(1X,2A,ES10.3)') Text,', max. abs. difference',Dev
Done = .true.
if (Dev > Thr) then
  write(u6,'(1X,3A,ES10.3)') 'Verification failed: ',Text,', max. abs. difference',Dev
  call Abend()
end if
end subroutine Report_Dev

!> Index and phase in the basis of the CI back end of the CSF iFile of the
!> split-graph order of the files (e.g. the reference CSFs of CIRO).
subroutine CI_File_Index_To_Basis(iState,LSM,n,iFile,iBasis,iPhase)
integer(kind=iwp), intent(in) :: iState, LSM, n, iFile
integer(kind=iwp), intent(out) :: iBasis, iPhase
real(kind=wp), allocatable :: E(:), F(:)
call mma_allocate(E,n,Label='E_File')
call mma_allocate(F,n,Label='F_Basis')
E(:) = Zero
E(iFile) = One
call CI_File_To_Basis(iState,LSM,n,E,F)
iBasis = maxloc(abs(F),1)
iPhase = nint(F(iBasis))
call mma_deallocate(F)
call mma_deallocate(E)
end subroutine CI_File_Index_To_Basis

!> GAS: the CSF map of symmetry LSM of state iState into its relaxed space.
subroutine SG_Rel_Map(iState,LSM)
use sguga, only: SG_Map_CSFs
integer(kind=iwp), intent(in) :: iState, LSM
integer(kind=iwp) :: n
if (allocated(Rel_Map(LSM)%M)) return
n = CIS(iState)%NCSF(LSM)
call mma_allocate(Rel_Map(LSM)%M,max(1,n),Label='Rel_Map')
Rel_Map(LSM)%M(:) = 0
call SG_Map_CSFs(iState,CI_St(iState)%iRel,LSM,n,Rel_Map(LSM)%M)
if (any(Rel_Map(LSM)%M(1:n) == 0)) then
  write(u6,*) 'SG_Rel_Map: a CSF of the GAS space is missing in the relaxed space'
  call Abend()
end if
end subroutine SG_Rel_Map

subroutine SG_Rel_Free()
integer(kind=iwp) :: k
do k=1,size(Rel_Map)
  call mma_deallocate(Rel_Map(k)%M,safe='*')
end do
end subroutine SG_Rel_Free

!> GAS: X (n CSFs of symmetry LSM of iState) into the relaxed space (XR,
!> zero elsewhere), and back (the components of the GAS CSFs).
subroutine SG_Rel_Embed(iState,LSM,n,X,XR)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: X(n)
real(kind=wp), intent(out) :: XR(:)
call SG_Rel_Map(iState,LSM)
XR(:) = Zero
XR(Rel_Map(LSM)%M(1:n)) = X(1:n)
end subroutine SG_Rel_Embed

subroutine SG_Rel_Project(iState,LSM,n,XR,X)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: XR(:)
real(kind=wp), intent(out) :: X(n)
call SG_Rel_Map(iState,LSM)
X(1:n) = XR(Rel_Map(LSM)%M(1:n))
end subroutine SG_Rel_Project

!> SGUGA sigma vector of state iState: sg_h_psi_ras (RAS, or Use_RAS) or
!> sg_h_psi (CAS); for GAS sg_h_psi_ras in the relaxed space.
subroutine SG_HPsi(iState,Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU,Use_RAS)
integer(kind=iwp), intent(in) :: iState, n, LSM, nTUVX, nTU
real(kind=wp), intent(in) :: Psi(n), TUVX(nTUVX), TU(nTU)
real(kind=wp), intent(out) :: Sigma(n)
logical(kind=iwp), intent(in) :: Use_RAS
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: PR(:), SR(:)
iR = CI_St(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(PR,nR,Label='PR')
  call mma_allocate(SR,nR,Label='SR')
  call SG_Rel_Embed(iState,LSM,n,Psi,PR)
  call sg_h_psi_ras(SGS(iR),CIS(iR),EXS(iR),PR,nR,LSM,SR,TUVX,nTUVX,TU,nTU)
  call SG_Rel_Project(iState,LSM,n,SR,Sigma)
  call mma_deallocate(SR)
  call mma_deallocate(PR)
else if (Use_RAS) then
  call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU)
else
  call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU)
end if
end subroutine SG_HPsi

!> SGUGA 2-RDM (sg_two_pdm) of state iState; for GAS in the relaxed space.
subroutine SG_Two_PDM_G(iState,C,n,LSM,P,PA,nP)
integer(kind=iwp), intent(in) :: iState, n, LSM, nP
real(kind=wp), intent(in) :: C(n)
real(kind=wp), intent(out) :: P(nP), PA(nP)
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: CR(:)
iR = CI_St(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(CR,nR,Label='CR')
  call SG_Rel_Embed(iState,LSM,n,C,CR)
  call sg_two_pdm(SGS(iR),CIS(iR),EXS(iR),CR,nR,LSM,P,PA,nP)
  call mma_deallocate(CR)
else
  call sg_two_pdm(SGS(iState),CIS(iState),EXS(iState),C,n,LSM,P,PA,nP)
end if
end subroutine SG_Two_PDM_G

!> SGUGA spin density (sg_spin_pdm, D the 1-RDM) of state iState; for GAS
!> in the relaxed space.
subroutine SG_Spin_PDM_G(iState,C,n,LSM,D,SD,nD)
integer(kind=iwp), intent(in) :: iState, n, LSM, nD
real(kind=wp), intent(in) :: C(n), D(nD)
real(kind=wp), intent(out) :: SD(nD)
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: CR(:)
iR = CI_St(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(CR,nR,Label='CR')
  call SG_Rel_Embed(iState,LSM,n,C,CR)
  call sg_spin_pdm(SGS(iR),CIS(iR),EXS(iR),CR,nR,LSM,D,SD,nD)
  call mma_deallocate(CR)
else
  call sg_spin_pdm(SGS(iState),CIS(iState),EXS(iState),C,n,LSM,D,SD,nD)
end if
end subroutine SG_Spin_PDM_G

!> SGUGA transition 1-RDM and spin 1-TDM (sg_trans_pdm) of state iState;
!> for GAS in the relaxed space.
subroutine SG_Trans_PDM_G(iState,B,K,n,LSM,nA,TS,TSS,Spin)
integer(kind=iwp), intent(in) :: iState, n, LSM, nA
real(kind=wp), intent(in) :: B(n), K(n)
real(kind=wp), intent(out) :: TS(nA*nA), TSS(nA*nA)
logical(kind=iwp), intent(in) :: Spin
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: BR(:), KR(:)
iR = CI_St(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(BR,nR,Label='BR')
  call mma_allocate(KR,nR,Label='KR')
  call SG_Rel_Embed(iState,LSM,n,B,BR)
  call SG_Rel_Embed(iState,LSM,n,K,KR)
  call sg_trans_pdm(SGS(iR),CIS(iR),EXS(iR),BR,KR,nR,LSM,nA,TS,TSS,Spin)
  call mma_deallocate(KR)
  call mma_deallocate(BR)
else
  call sg_trans_pdm(SGS(iState),CIS(iState),EXS(iState),B,K,n,LSM,nA,TS,TSS,Spin)
end if
end subroutine SG_Trans_PDM_G

!> Exact SGUGA diagonal <I|H|I> of state iState (Mode 'FAST': sg_hdiag_fast,
!> 'REF': sg_hdiag_ref); for GAS in the relaxed space.
subroutine SG_HDiag_G(Mode,iState,n,LSM,TUVX,nTUVX,TU,nTU,Diag)
character(len=*), intent(in) :: Mode
integer(kind=iwp), intent(in) :: iState, n, LSM, nTUVX, nTU
real(kind=wp), intent(in) :: TUVX(nTUVX), TU(nTU)
real(kind=wp), intent(out) :: Diag(n)
integer(kind=iwp) :: iR, iS, nR
real(kind=wp), allocatable :: DR(:)
iR = CI_St(iState)%iRel
iS = iState
nR = n
if (iR > 0) then
  iS = iR
  nR = CIS(iR)%NCSF(LSM)
end if
call mma_allocate(DR,nR,Label='DR')
if (Mode == 'REF') then
  call sg_hdiag_ref(SGS(iS),CIS(iS),EXS(iS),nR,LSM,TUVX,nTUVX,TU,nTU,DR)
else
  call sg_hdiag_fast(SGS(iS),CIS(iS),EXS(iS),nR,LSM,TUVX,nTUVX,TU,nTU,DR)
end if
if (iR > 0) then
  call SG_Rel_Project(iState,LSM,n,DR,Diag)
else
  Diag(:) = DR(:)
end if
call mma_deallocate(DR)
end subroutine SG_HDiag_G

!> SGUGA covers the CI space of RASSCF: CAS and RAS directly, GAS through
!> the relaxed space (set up by SG_Setup_RASSCF; HEXS/DEXS are masks).
function SG_Avail() result(Avail)
integer(kind=iwp), parameter :: iState = 1
logical(kind=iwp) :: Avail
Avail = (.not. iDoGAS) .or. (CI_St(iState)%iRel > 0)
end function SG_Avail

!> HEXS/DEXS: is a configuration with nOcc(k) electrons in partition k of
!> Space_CI removed (the removal rules of the CI space)?
function HEXS_Removed(nOcc) result(Rem)
use CI_Space_Def, only: Space_Allowed
integer(kind=iwp), intent(in) :: nOcc(:)
logical(kind=iwp) :: Rem
Rem = .not. Space_Allowed(Space_CI,nOcc(1:Space_CI%nPart))
end function HEXS_Removed

!> HEXS/DEXS: the kept CSFs of symmetry STSYM in the split-graph basis,
!> from the SGUGA configurations.
subroutine HEXS_Mask_SG(iState,Keep)
use sguga, only: SG_Mk_Cnf_List, SG_Cnf_Occ
integer(kind=iwp), intent(in) :: iState
logical(kind=iwp), intent(out) :: Keep(:)
type(SGCnfList) :: Cnf
integer(kind=iwp) :: k, iLev, nOcc(MxGAS)
integer(kind=iwp), allocatable :: Occ(:)
call SG_Mk_Cnf_List(iState,CI_St(iState)%STSym,Cnf)
call mma_allocate(Occ,max(1,Cnf%nLev),Label='HEXS_Occ')
Keep(:) = .true.
do k=1,Cnf%nCnf
  call SG_Cnf_Occ(Cnf,k,Occ)
  nOcc(:) = 0
  do iLev=1,Cnf%nLev
    nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+Occ(iLev)
  end do
  if (HEXS_Removed(nOcc)) Keep(Cnf%iCSF(Cnf%iOff(k)+1:Cnf%iOff(k+1))) = .false.
end do
call mma_deallocate(Occ)
call SG_Free_Cnf_List(Cnf)
end subroutine HEXS_Mask_SG

!> HEXS/DEXS: the kept CSFs in the basis of Faroald (citrans), from its
!> configurations (the orbitals of Faroald are the levels).
subroutine HEXS_Mask_Faro(Keep)
logical(kind=iwp), intent(out) :: Keep(:)
integer(kind=iwp) :: k, k1, nc, nd, iLev, nOcc(MxGAS)
integer(kind=iwp), allocatable :: da(:), db(:)
real(kind=wp), allocatable :: cf(:,:)
Keep(:) = .true.
k = 1
do while (k <= nConf)
  call citrans_conf_dets(k,k1,nc,nd,da,db,cf)
  nOcc(:) = 0
  do iLev=1,NAC
    if (btest(da(1),iLev-1)) nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+1
    if (btest(db(1),iLev-1)) nOcc(Space_CI%LevPart(iLev)) = nOcc(Space_CI%LevPart(iLev))+1
  end do
  if (HEXS_Removed(nOcc)) Keep(k1:k1+nc-1) = .false.
  k = k1+nc
  call mma_deallocate(da)
  call mma_deallocate(db)
  call mma_deallocate(cf)
end do
end subroutine HEXS_Mask_Faro

!> HEXS/DEXS: the masks of the kept CSFs (end of CI_Initialize). Without
!> Lucia the number of kept CSFs (CI_St%nCSF_Kept) and the cap of NSEL come from
!> here; with Lucia they are compared.
subroutine HEXS_Setup(iState,Do_SGUGA)
use rasscf_global, only: NSEL
integer(kind=iwp), intent(in) :: iState
logical(kind=iwp), intent(in) :: Do_SGUGA
integer(kind=iwp) :: nKeep
call mma_deallocate(HEXS_Keep,safe='*')
call mma_deallocate(HEXS_Keep_SG,safe='*')
! the removals are part of the CI space (CI_Space_Setup)
if (Space_CI%nRem <= 0) return
if (Do_SGUGA .and. SG_Avail()) then
  call mma_allocate(HEXS_Keep_SG,nConf,Label='HEXS_Keep_SG')
  call HEXS_Mask_SG(iState,HEXS_Keep_SG)
end if
if (SG_CI()) then
  call mma_allocate(HEXS_Keep,nConf,Label='HEXS_Keep')
  HEXS_Keep(:) = HEXS_Keep_SG(:)
else if (DoFaro) then
  call mma_allocate(HEXS_Keep,nConf,Label='HEXS_Keep')
  call HEXS_Mask_Faro(HEXS_Keep)
end if
if (allocated(HEXS_Keep)) then
  nKeep = count(HEXS_Keep)
else if (allocated(HEXS_Keep_SG)) then
  nKeep = count(HEXS_Keep_SG)
else
  return
end if
write(u6,'(1X,A,I0,A,I0,A)') 'HEXS/DEXS: ',nKeep,' of ',nConf,' CSFs kept'
if (.not. Lucia_Active) then
  CI_St(iState)%nCSF_Kept = nKeep
  if (NSEL > nKeep) NSEL = nKeep
else if (nKeep /= CI_St(iState)%nCSF_Kept) then
  write(u6,'(1X,A,I0)') 'HEXS/DEXS: Lucia keeps (nCSF_HEXS) ',CI_St(iState)%nCSF_Kept
end if
end subroutine HEXS_Setup

!> HEXS/DEXS: v := P v (zero on the removed CSFs), if the mask is set and
!> matches v.
subroutine HEXS_Proj(Keep,v)
logical(kind=iwp), allocatable, intent(in) :: Keep(:)
real(kind=wp), intent(inout) :: v(:)
if (.not. allocated(Keep)) return
if (size(Keep) /= size(v)) return
where (.not. Keep) v = Zero
end subroutine HEXS_Proj

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Sigma vector Sigma_Vec = H CI_Vec of the CI back end; with HEXS/DEXS
!> P H P (Mk_H_Psi_P on the projected vector, the result projected).
Subroutine Mk_H_Psi(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,nTU,TU,nTUVX,TUVX)
integer(kind=iwp), intent(in):: iState, STSYM, nCSF, ntemp, nTU, nTUVX
real(kind=wp), intent(in) :: CI_Vec(nCSF), TU(nTU), TUVX(nTUVX)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
real(kind=wp), allocatable :: CM(:)
if (.not. allocated(HEXS_Keep)) then
  call Mk_H_Psi_P(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,nTU,TU,nTUVX,TUVX)
  return
end if
call mma_allocate(CM,nCSF,Label='HEXS_CM')
CM(:) = CI_Vec(:)
call HEXS_Proj(HEXS_Keep,CM)
call Mk_H_Psi_P(iState,STSYM,nCSF,CM,Sigma_Vec,ctemp,sigtemp,ntemp,nTU,TU,nTUVX,TUVX)
call HEXS_Proj(HEXS_Keep,Sigma_Vec)
call mma_deallocate(CM)
end subroutine Mk_H_Psi

Subroutine Mk_H_Psi_P(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,nTU,TU,nTUVX,TUVX)
use lucia_data, only: Sigma_on_disk, IREFSM
use cands, only: ICSM, ISSM
Implicit None

integer(kind=iwp), intent(in):: iState, STSYM, nCSF
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
integer(kind=iwp), intent(in) :: ntemp
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

integer(kind=iwp) :: IREFSM_SAVE, nDet

real(kind=wp) :: Check_Href, Check_H, DMax
logical(kind=iwp) :: Lucia_Done, SG_Done
real(kind=wp), allocatable ::  SG_PSI(:), SG_SGM(:), SG_RAS(:)

! SGCI: the CI is in the SGUGA (split-graph) CSF basis. Under
! MOLCAS_CI_VERIFY=SIGMA the sigma vector is compared with that of Lucia
! (vector reordered to the Lucia basis and the result back).
if (SG_CI()) then
  call SG_HPsi(iState,CI_Vec,nCSF,STSYM,Sigma_Vec,TUVX,nTUVX,TU,nTU,CI_St(iState)%nRsPrt /= 1)
  call HEXS_Proj(HEXS_Keep,Sigma_Vec)
  if (CI_Verify('SIGMA')) then
    call mma_allocate(SG_PSI,max(nCSF,CI_St(iState)%nDet(STSYM)),Label='SG_PSI')
    call mma_allocate(SG_SGM,max(nCSF,CI_St(iState)%nDet(STSYM)),Label='SG_SGM')
    call SG_Reord(iState,STSYM,1,nCSF,CI_Vec,SG_PSI)
    call Lucia_Sigma(SG_PSI,SG_SGM)
    call SG_Reord(iState,STSYM,0,nCSF,SG_SGM,SG_PSI)
    DMax = maxval(abs(SG_PSI(1:nCSF)-Sigma_Vec(:)))
    if (.not. HPsi_SG_Verified) then
      write(u6,'(1X,A,ES10.3)') 'Mk_H_Psi (SGCI): SGUGA sigma vs Lucia, max. abs. difference',DMax
      HPsi_SG_Verified = .true.
    end if
    if (DMax > 1.0e-10_wp*max(One,maxval(abs(Sigma_Vec(:))))) then
      write(u6,'(1X,A,ES10.3)') 'SGUGA error in H|Psi> (SGCI), max. abs. difference',DMax
      call Abend()
    end if
    call mma_deallocate(SG_SGM)
    call mma_deallocate(SG_PSI)
  end if
  return
end if

if (DOFARO) then

  call Faroald_Ints(TU,TUVX)
  call Faroald_Sigma(nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb)
  call HEXS_Proj(HEXS_Keep,Sigma_Vec)
  ! CI development (MOLCAS_CI_VERIFY=SIGMA): the Faroald sigma vector
  ! against Lucia (and below against SGUGA)
  if (CI_Verify('SIGMA')) then
    call mma_allocate(SG_SGM,nCSF,Label='Lucia_SGM')
    call Ref_Sigma('Lucia',iState,STSYM,nCSF,CI_Vec,SG_SGM,nTU,TU,nTUVX,TUVX)
    DMax = maxval(abs(SG_SGM(:)-Sigma_Vec(:)))
    Lucia_Done = HPsi_Faro_Verified
    call Report_Dev('Mk_H_Psi (FARO): Faroald sigma vs Lucia',DMax,1.0e-10_wp*max(One,maxval(abs(Sigma_Vec(:)))),Lucia_Done)
    call mma_deallocate(SG_SGM)
  end if

else

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  nDET = CI_St(iState)%nDet(IREFSM)

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ctemp(1:nCSF) = CI_Vec(1:nCSF)
  sigtemp(:) = Zero
  call csdtvc(ctemp,sigtemp,1,STSym,0)
  ctemp(1:nDet)=sigtemp(1:nDet)

  ! Calling Lucia to determine the sigma vector
  call Lucia_Util('Sigma',                   &
                  CI_Vector=ctemp(:),        &
                  Sigma_Vector=sigtemp(:),   &
                  nTU=Size(TU),TU=TU,        &
                  nTUVX=Size(TUVX),TUVX=TUVX)

  ! Set mark so densi_master knows that the Sigma-vector exists on disk.
  Sigma_on_disk = .true.
  ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
  ! stored in Sigma_vec.
  call CSDTVC(Sigma_Vec,sigtemp,2,STSym,0)
  sigtemp(1:nCSF)=Sigma_Vec(1:nCSF)

  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM

end if

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the SGUGA sigma
! vector, which is then used. With RAS partitions it comes from
! sg_h_psi_ras (no intermediate states outside the RAS space); for CAS from
! sg_h_psi, checked once per run against sg_h_psi_ras. With HEXS/DEXS
! projected (P H P, HEXS_Keep_SG).
If (CI_Verify('SIGMA') .and. SG_Avail()) Then
   Call mma_allocate(SG_PSI,nCSF,Label='SG_PSI')
   call CI_Basis_To_File(iState,STSYM,nCSF,CI_VEC,SG_PSI)
   call HEXS_Proj(HEXS_Keep_SG,SG_PSI)
   Check_Href=CheckSum(Sigma_Vec,nCSF)
   Call mma_allocate(SG_SGM,nCSF,Label='SG_SGM')
   if ((CI_St(iState)%nRsPrt == 1) .and. (.not. iDoGAS)) then
     Call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_SGM,TUVX,nTUVX,TU,nTU)
     if (.not. HPsi_RAS_Verified) then
       call mma_allocate(SG_RAS,nCSF,Label='SG_RAS')
       call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_RAS,TUVX,nTUVX,TU,nTU)
       DMax = maxval(abs(SG_RAS(:)-SG_SGM(:)))
       write(u6,'(1X,A,ES10.3)') 'Mk_H_Psi: sg_h_psi_ras vs sg_h_psi, max. abs. difference',DMax
       if (DMax > 1.0e-12_wp*max(One,maxval(abs(SG_SGM(:))))) then
         write(u6,*) 'SGUGA error in sg_h_psi_ras'
         call Abend()
       end if
       call mma_deallocate(SG_RAS)
       HPsi_RAS_Verified = .true.
     end if
   else
     call SG_HPsi(iState,SG_PSI,nCSF,STSYM,SG_SGM,TUVX,nTUVX,TU,nTU,.true.)
   end if
   call HEXS_Proj(HEXS_Keep_SG,SG_SGM)
   call CI_File_To_Basis(iState,STSYM,nCSF,SG_SGM,SG_PSI)
   if ((CI_St(iState)%nRsPrt /= 1) .or. iDoGAS) then
     DMax = maxval(abs(SG_PSI(:)-Sigma_Vec(:)))
     if (.not. HPsi_RAS_Verified) &
       write(u6,'(1X,A,ES10.3)') 'Mk_H_Psi: sg_h_psi_ras vs CI back end (RAS), max. abs. difference',DMax
     HPsi_RAS_Verified = .true.
     if (DMax > 1.0e-10_wp*max(One,maxval(abs(Sigma_Vec(:))))) then
       write(u6,'(1X,A,ES10.3)') 'SGUGA error in sg_h_psi_ras (RAS), max. abs. difference',DMax
       call Abend()
     end if
   end if
   if (DoFaro) then
     SG_Done = HPsi_Faro_Verified
     call Report_Dev('Mk_H_Psi (FARO): Faroald sigma vs SGUGA',maxval(abs(SG_PSI(:)-Sigma_Vec(:))), &
                     1.0e-10_wp*max(One,maxval(abs(Sigma_Vec(:)))),SG_Done)
     HPsi_Faro_Verified = .true.
   end if
   Check_H   =CheckSum(SG_PSI,nCSF)
   If (Abs(Check_HRef-Check_H)/nCSF>1.0E-12_wp) Then
      Write (u6,*) 'SGUGA error in H|Psi>'
      Write (u6,*) 'Check_HRef=',Check_HRef
      Write (u6,*) 'Check_H   =',Check_H
      Call RecPrt('Sigma_Vec',' ',Sigma_Vec,1,nCSF)
      Call RecPrt('SG_SGM',' ',SG_PSI,1,nCSF)
      Call Abend()
   End If
   Sigma_Vec(1:nCSF) = SG_PSI(1:nCSF)
   Call mma_deallocate(SG_SGM)
   Call mma_deallocate(SG_PSI)
End if

contains

  ! Lucia sigma vector of the CSF vector Cin (Lucia basis), the reference
  ! of SGCI; Sout must hold max(nCSF,nDet) elements (CSDTVC).
  subroutine Lucia_Sigma(Cin,Sout)
    real(kind=wp), intent(in) :: Cin(*)
    real(kind=wp), intent(inout) :: Sout(*)
    integer(kind=iwp) :: IREFSM_Keep, nD

    IREFSM_Keep = IREFSM
    IREFSM = STSYM
    ICSM = STSYM
    ISSM = STSYM
    nD = CI_St(iState)%nDet(IREFSM)
    ctemp(1:nCSF) = Cin(1:nCSF)
    sigtemp(:) = Zero
    call csdtvc(ctemp,sigtemp,1,STSym,0)
    ctemp(1:nD) = sigtemp(1:nD)
    call Lucia_Util('Sigma',CI_Vector=ctemp(:),Sigma_Vector=sigtemp(:),nTU=Size(TU),TU=TU,nTUVX=Size(TUVX),TUVX=TUVX)
    Sigma_on_disk = .true.
    call CSDTVC(Sout,sigtemp,2,STSym,0)
    IREFSM = IREFSM_Keep
    ICSM = IREFSM
    ISSM = IREFSM

  end subroutine Lucia_Sigma

End Subroutine Mk_H_Psi_P

!> Unpack the active one- and two-electron integrals (triangular packing)
!> into the square arrays htu and gtuvx of Faroald.
!> The integrals of Faroald (htu, gtuvx over its orbitals = the levels)
!> from TU/TUVX (triangular, symmetry-blocked active orbitals), via Faro_L2O.
subroutine Faroald_Ints(TU,TUVX)
real(kind=wp), intent(in) :: TU(:), TUVX(:)
integer(kind=iwp) :: it, iu, iv, ix, ot, ou, ov, ox, otu, ovx

htu(:,:) = Zero
gtuvx(:,:,:,:) = Zero
do it=1,my_norb
  ot = Faro_L2O(it)
  do iu=1,my_norb
    ou = Faro_L2O(iu)
    otu = iTri_F(ot,ou)
    htu(it,iu) = TU(otu)
    do iv=1,my_norb
      ov = Faro_L2O(iv)
      do ix=1,my_norb
        ox = Faro_L2O(ix)
        ovx = iTri_F(ov,ox)
        gtuvx(it,iu,iv,ix) = TUVX(iTri_F(otu,ovx))
      end do
    end do
  end do
end do

contains

  pure function iTri_F(i,j) result(ij)
    integer(kind=iwp), intent(in) :: i, j
    integer(kind=iwp) :: ij
    ij = max(i,j)*(max(i,j)-1)/2+min(i,j)
  end function iTri_F

end subroutine Faroald_Ints

!> Faroald sigma vector Sigma_Vec = H CI_Vec in its own CSF basis (citrans), with the
!> integrals set by Faroald_Ints. ctemp and sigtemp are scratch of at least
!> nDetA*nDetB elements (Faroald determinant basis).
subroutine Faroald_Sigma(nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb)
#ifdef _Compute_Full_CI_Matrix
use SGUGA_Sigma, only: sg_ci_full
#endif
integer(kind=iwp), intent(in) :: nCSF, ntemp, ndeta, ndetb
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
real(kind=wp), pointer :: Faroald_PSI(:,:), Faroald_SGM(:,:)

real(kind=wp), allocatable :: VRas(:), SRas(:), SFull(:)
real(kind=wp), allocatable, target :: PFull(:), QFull(:)
real(kind=wp) :: DMax

! RAS: string-driven, only the determinants of the RAS space
if (ras_active) then
  call mma_allocate(VRas,max(1,n_ras),Label='VRas')
  call mma_allocate(SRas,max(1,n_ras),Label='SRas')
  call CITRANS_CSF2RAS(CI_Vec,VRas)
  call RAS_SIGMA(HTU,GTUVX,VRas,SRas)
  call CITRANS_RAS2CSF(SRas,Sigma_Vec)
  call mma_deallocate(VRas)
  call mma_deallocate(SRas)
  ! CI development (MOLCAS_CI_VERIFY=SIGMA): once against the sigma of the
  ! full determinant space (P H P), if its tables are set up (faroald_full)
  if (HPsi_FaroRAS_Verified .or. (.not. CI_Verify('SIGMA')) .or. (.not. faroald_full)) return
  call mma_allocate(SFull,nCSF,Label='SFull')
  call mma_allocate(PFull,nDetA*nDetB,Label='PFull')
  call mma_allocate(QFull,nDetA*nDetB,Label='QFull')
  Faroald_Psi(1:nDetA,1:nDetB) => PFull(:)
  Faroald_SGM(1:nDetA,1:nDetB) => QFull(:)
else
  Faroald_Psi(1:nDetA,1:nDetB) => ctemp(:)
  Faroald_SGM(1:nDetA,1:nDetB) => sigtemp(:)
end if

! (the CI vectors are in the configuration order of Faroald, citrans)
Faroald_PSI(:,:) = Zero
call CITRANS_CSF2SD(CI_Vec,Faroald_PSI)
Faroald_SGM(:,:) = Zero
call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)

Faroald_Psi => Null()
Faroald_SGM => Null()

if (ras_active) then
  call CITRANS_SD2CSF(QFull,SFull)
  DMax = maxval(abs(SFull(:)-Sigma_Vec(:)))
  call Report_Dev('Faroald_Sigma: string-driven RAS sigma vs full determinant space',DMax, &
                  1.0e-10_wp*max(One,maxval(abs(SFull(:)))),HPsi_FaroRAS_Verified)
  call mma_deallocate(SFull)
  call mma_deallocate(PFull)
  call mma_deallocate(QFull)
else
  call CITRANS_SD2CSF(sigtemp,Sigma_Vec)
end if

end subroutine Faroald_Sigma

Subroutine Mk_H_Psi_Det(istate,STSYM,nDet,CI_Vec,Sigma_Vec,nTU,TU,nTUVX,TUVX)
use lucia_data, only: Sigma_on_disk, IREFSM
use sd_csf_codec_module, only: sd_csf_codec, sd_to_csf, csf_to_sd
use cands, only: ICSM, ISSM
Implicit None

integer(kind=iwp), intent(in):: istate, STSYM, nDet
real(kind=wp), intent(inout),target :: CI_Vec(nDet)
real(kind=wp), intent(inout),target :: Sigma_Vec(nDet)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

real(kind=wp), allocatable, target:: temp1(:), temp2(:)  ! Local scratch
real(kind=wp), pointer:: Faroald_PSI(:,:), Faroald_SGM(:,:)

integer(kind=iwp) :: nCSF
integer(kind=iwp) :: iDet
integer(kind=iwp) :: IREFSM_SAVE

#ifdef _Compute_Full_CI_Matrix
real(kind=wp), allocatable :: Full(:,:)
#endif

!Note that we Lucia code doesn't explicily spin adapt the SD. Rather it just process the input coefficients.
! If we project the starting coefficients to the spin-adapted iCSF basis and then expand it back to the SD basis
! we have effectivly generated a SD coefficient vector that is implicitly spin-adapted. Since the CASVB code generates
! starting CI vectors that are not spin-adapted there will be a difference in the convergence patterns.

#ifdef _Compute_Full_CI_Matrix
! Reference computation of the full CI matrix.
! in CSF basis
nCSF=CIS(istate)%nCSF(STSYM)
Call mma_allocate(Full,nCSF,nCSF,Label='Full')
Call sg_ci_full(SGS(istate),CIS(istate),EXS(istate),nCSF,STSym, &
                TUVX,nTUVX, &
                TU,nTU, &
                Full)
Call mma_deallocate(Full)

Call mma_allocate(Full,nDet,nDet,Label='Full')
Call ci_full_det(nDet,STSym, &
                TUVX,nTUVX, &
                TU,nTU, &
                Full)
Call mma_deallocate(Full)
! in SD basis
#endif

if (DOFARO) then

  call Faroald_Ints(TU,TUVX)

  ! We need to convert this to the correct form! The input is in LUCIA format.

  nCSF=CIS(istate)%nCSF(STSym) ! ToDo: This depends on the SGUGA part of the code.
  if (nCSF>nDet) Then
     Write(u6,*) 'Mk_H_Psi_Det: nCSF>nDet'
     Call Abend()
  End If

  Call mma_allocate(temp1,nDet,Label='temp1')
  Call mma_allocate(temp2,nDet,Label='temp2')
  Faroald_Psi(1:nDetA,1:nDetB) => temp2(:)
  Faroald_SGM(1:nDetA,1:nDetB) => temp1(:)

  temp2(:)=CI_Vec(:)
  ! Convert from LUCIA SD basis to Lucia CSF basis
  temp1(:)=Zero
  Call SD_CSF_Codec(temp2,temp1,SD_to_CSF,STSym)
  ! Now use the FAROALD machinary to go to FAROALD SD basis
  call SG_REORD(iState,STSYM,0,nCSF,temp1,temp2)
  call CITRANS_SORT('C',temp2,temp1)
  Faroald_PSI(:,:) = Zero
  call CITRANS_CSF2SD(temp1,Faroald_PSI)


  call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)

  call CITRANS_SD2CSF(Faroald_SGM,Sigma_Vec)
  call CITRANS_SORT('O',Sigma_Vec,temp2)
  call SG_Reord(iState,STSYM,1,nCSF,temp2,temp1)

  ! Convert from LUCIA CSF basis to LUCIA SD basis
  Call SD_CSF_Codec(Sigma_Vec,temp1,CSF_to_SD,STSym)

  Faroald_Psi => Null()
  Faroald_SGM => Null()
  Call mma_deallocate(temp2)
  Call mma_deallocate(temp1)

else

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ! Calling Lucia to determine the sigma vector
  call Lucia_Util('Sigma',                   &
                  CI_Vector=CI_Vec(:),        &
                  Sigma_Vector=Sigma_Vec(:),   &
                  nTU=Size(TU),TU=TU,        &
                  nTUVX=Size(TUVX),TUVX=TUVX)

  ! Set mark so densi_master knows that the Sigma-vector exists on disk.
  Sigma_on_disk = .true.
  ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
  ! stored in Sigma_vec.

  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM
endif

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the SGUGA sigma
! vector (CASSCF without symmetry only), which is then used. The
! determinant (CASVB) route needs the SD/CSF codec, which exists for CAS
! only; with RAS partitions the check is skipped.
If (CI_Verify('SIGMA')) Then
  If (SGS(iState)%nRsPrt /= 1) Then
     Return
  End If
  If (SGS(iState)%nSym /= 1) Return
  nCSF=CIS(istate)%nCSF(STSym)
  Call mma_allocate(temp1,nDet,Label='temp1')
  Call mma_allocate(temp2,nDet,Label='temp2')

  temp2(:)=CI_Vec(:)
  ! Convert from LUCIA SD basis to Lucia CSF basis
  temp1(:)=Zero
  Call SD_CSF_Codec(temp2,temp1,SD_to_CSF,STSYM)
  ! Convert from Lucia CSF basis to SGUGA CSF basis
  !     iMode=0 --> from SGA to split graph GUGA order                   *
  !     iMode=1 --> from split graph GUGA to SGA order                   *
  call SG_Reord(iState,STSYM,0,nCSF,temp1,temp2)

  call SG_HPsi(iState,temp2,nCSF,STSYM,temp1,TUVX,nTUVX,TU,nTU,CI_St(iState)%nRsPrt /= 1)

  ! Convert from SGUGA CSF basis to Lucia CSF basis
  call SG_Reord(iState,STSYM,1,nCSF,temp1,temp2)
  ! Convert from Lucia CSF basis to Lucia SD Basis
   call SD_CSF_Codec(temp1,temp2,CSF_to_SD,STSym)

  Do iDet = 1, nDet
     If (Abs(temp1(iDet)-Sigma_vec(iDet))>1.0e-10_wp) Then
        Write (u6,*) 'Mk_H_Psi_Det: Abs(temp1(iDet)-Sigma_vec(iDet))>1.0e-10_wp)'
        Write (u6,*) ' iDet=',iDet
        Write (u6,*) 'temp1(iDet)=    ',temp1(iDet)
        Write (u6,*) 'Sigma_vec(iDet)=', Sigma_vec(iDet)
        Call Abend()
     End If
  End Do

  Sigma_Vec(1:nDet)=temp1(1:nDet)
  Call mma_deallocate(temp1)
  Call mma_deallocate(temp2)
End If

#ifdef _Compute_Full_CI_Matrix
Contains
subroutine ci_full_det(nDet,STSym, &
                       TUVX,nTUVX, &
                       TU,nTU, &
                       Full)
use constants, only: Zero, One
use lucia_data, only: IREFSM
use cands, only: ICSM, ISSM

  integer(kind=iwp), intent(in):: nDet,STSym,nTUVX,nTU
  real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout):: Full(nDet,nDet)

  integer(kind=iwp) IREFSM_SAVE, iDet, jDet

  real(kind=wp), allocatable :: CI_VEC(:), Sigma_VEC(:)

  Call mma_allocate(CI_Vec,nDet,Label='CI_Vec')
  Call mma_allocate(Sigma_Vec,nDet,Label='Sigma_Vec')

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  Do iDet = 1, nDet
     CI_Vec(:) = Zero
     CI_Vec(iDet) = One

  ! Calling Lucia to determine the sigma vector
     Sigma_Vec(:) = Zero
     call Lucia_Util('Sigma',                   &
                     CI_Vector=CI_Vec(:),        &
                     Sigma_Vector=Sigma_Vec(:),   &
                     nTU=nTU,TU=TU,        &
                     nTUVX=nTUVX,TUVX=TUVX)

     Do jDet = 1, nDet
        Full(iDet,jDet) = Sigma_Vec(jDet)
     End Do

  End Do
  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM

  Call mma_deallocate(CI_Vec)
  Call mma_deallocate(Sigma_Vec)

  Call RecPrt('CI Matrix in SD basis',' ',Full,ndet,ndet)

end subroutine ci_full_det
#endif

End Subroutine Mk_H_Psi_Det

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_T1DM(Bra_Vec, Ket_Vec, nVec, T1DM, nT1DM, T1SDM)
use Lucia_Data, only: DTmp, DSTmp
implicit none
integer(kind=iwp), intent(in) :: nVec, nT1DM
real(kind=wp), intent(_IN_):: Bra_Vec(nVec), Ket_Vec(nVec)
real(kind=wp), intent(out):: T1DM(nT1DM)
real(kind=wp), optional, intent(out):: T1SDM(nT1DM)

real(kind=wp), allocatable:: T1SDM_loc(:), Bra_SD(:,:), Ket_SD(:,:)
integer(kind=iwp), parameter :: iState=1
logical(kind=iwp) :: free_D, free_DS

! SGCI: transition densities from SGUGA in its basis; under SIGMA compared
! with Lucia (both vectors reordered to the Lucia basis).
if (SG_CI()) then
  call SGCI_T1DM()
  return
end if

If (DoFaro) Then
   Call mma_allocate(T1SDM_loc,nT1DM,Label='T1SDM')

   If (ras_active) Then
     ! RAS: in the RAS string classes (faroald_ras)
     Call mma_allocate(Bra_SD,max(1,n_ras),1,Label='Psi')
     Call mma_allocate(Ket_SD,max(1,n_ras),1,Label='Psi')
     call CITRANS_CSF2RAS(Bra_Vec,Bra_SD(:,1))
     call CITRANS_CSF2RAS(Ket_Vec,Ket_SD(:,1))
     Call RAS_Trans_One_PDM(Ket_SD(:,1),Bra_SD(:,1),T1DM,T1SDM_loc)
   Else
     Call mma_allocate(Bra_SD,nDetA,nDetB,Label='Psi')
     Call mma_allocate(Ket_SD,nDetA,nDetB,Label='Psi')

     ! (the CI vectors are in the configuration order of Faroald, citrans)
     Bra_SD(:,:)=Zero
     call CITRANS_CSF2SD(Bra_Vec,Bra_SD)
     Ket_SD(:,:)=Zero
     call CITRANS_CSF2SD(Ket_Vec,Ket_SD)

!    Don't ask me why we need to do it the other way?
!    Call Transition_One_PDM(Bra_SD,Ket_SD,T1DM,T1SDM_loc)
     Call Transition_One_PDM(Ket_SD,Bra_SD,T1DM,T1SDM_loc)
   End If

   If (present(T1SDM)) T1SDM(1:nT1DM) = T1SDM_loc(1:nT1DM)
   If (CI_Verify('SIGMA')) call Faroald_Check_T1DM()

   Call mma_deallocate(Ket_SD)
   Call mma_deallocate(Bra_SD)
   Call mma_deallocate(T1SDM_loc)
Else
   Free_D =.False.
   Free_DS=.False.
   If (.Not.Allocated(DTMP)) Then
      Free_D =.True.
      Call mma_allocate(DTMP,nT1DM,Label='DTMP')
   End If
   If (.Not.Allocated(DSTMP)) Then
      Free_DS=.True.
      Call mma_allocate(DSTMP,nT1DM,Label='DSTMP')
   End If

   call Lucia_Util('Densi',CI_Vector=Ket_Vec(:),RVec=Bra_Vec(:))
   T1DM(1:nT1DM)=DTMP(1:nT1DM)
   If (present(T1SDM)) T1SDM(1:nT1DM) = DSTMP(1:nT1DM)

   If (Free_D ) Call mma_deallocate(DTMP )
   If (Free_DS) Call mma_deallocate(DSTMP)
End If

contains

  ! FARO (CI development, SIGMA): the Faroald transition densities against
  ! Lucia (no spin density for singlets) and SGUGA (sg_trans_pdm, whose
  ! layout is the transpose).
  subroutine Faroald_Check_T1DM()
    real(kind=wp), allocatable :: BL(:), KL(:), TS(:), TSS(:)
    real(kind=wp) :: Dev(2)
    integer(kind=iwp) :: p, q
    logical(kind=iwp) :: Done
    real(kind=wp), parameter :: ThrD = 1.0e-10_wp

    call mma_allocate(BL,nVec,Label='BL')
    call mma_allocate(KL,nVec,Label='KL')
    call CI_Basis_To_Lucia(iState,CI_St(iState)%STSym,nVec,Bra_Vec,BL)
    call CI_Basis_To_Lucia(iState,CI_St(iState)%STSym,nVec,Ket_Vec,KL)
    call mma_allocate(DTMP,max(nT1DM,NAC*NAC),Label='DTMP')
    call mma_allocate(DSTMP,max(nT1DM,NAC*NAC),Label='DSTMP')
    call Lucia_Util('Densi',CI_Vector=KL(:),RVec=BL(:))
    Dev(1) = maxval(abs(DTMP(1:NAC*NAC)-T1DM(1:NAC*NAC)))
    Dev(2) = Zero
    if (CI_St(iState)%iSpin > 1) Dev(2) = maxval(abs(DSTMP(1:NAC*NAC)-T1SDM_loc(1:NAC*NAC)))
    Done = T1DM_Faro_Verified
    call Report_Dev('Mk_T1DM (FARO): Faroald T1DM, T1SDM vs Lucia',maxval(Dev),ThrD,Done)
    call mma_deallocate(DSTMP)
    call mma_deallocate(DTMP)
    if (SG_Avail()) then
      call CI_Basis_To_File(iState,CI_St(iState)%STSym,nVec,Bra_Vec,BL)
      call CI_Basis_To_File(iState,CI_St(iState)%STSym,nVec,Ket_Vec,KL)
      call mma_allocate(TS,NAC*NAC,Label='TS')
      call mma_allocate(TSS,NAC*NAC,Label='TSS')
      call SG_Trans_PDM_G(iState,BL,KL,nVec,CI_St(iState)%STSym,NAC,TS,TSS,.true.)
      Dev(:) = Zero
      do q=1,NAC
        do p=1,NAC
          Dev(1) = max(Dev(1),abs(T1DM(p+(q-1)*NAC)-TS(q+(p-1)*NAC)))
          Dev(2) = max(Dev(2),abs(T1SDM_loc(p+(q-1)*NAC)-TSS(q+(p-1)*NAC)))
        end do
      end do
      Done = T1DM_Faro_Verified
      call Report_Dev('Mk_T1DM (FARO): Faroald T1DM, T1SDM vs SGUGA',maxval(Dev),ThrD,Done)
      call mma_deallocate(TSS)
      call mma_deallocate(TS)
    end if
    T1DM_Faro_Verified = .true.
    call mma_deallocate(KL)
    call mma_deallocate(BL)
  end subroutine Faroald_Check_T1DM

  subroutine SGCI_T1DM()
    real(kind=wp), allocatable :: TS(:), TSS(:), BL(:), KL(:)
    real(kind=wp) :: Dev(2)
    integer(kind=iwp) :: p, q
    logical(kind=iwp) :: Spin
    real(kind=wp), parameter :: ThrD = 1.0e-10_wp

    Spin = present(T1SDM)
    call mma_allocate(TS,NAC*NAC,Label='TS')
    call mma_allocate(TSS,NAC*NAC,Label='TSS')
    call SG_Trans_PDM_G(iState,Bra_Vec,Ket_Vec,nVec,CI_St(iState)%STSym,NAC,TS,TSS,Spin)
    ! layout of the CI interface (as Lucia): element p+(q-1)*NAC is
    ! <Bra|E_qp|Ket>, the transpose of sg_trans_pdm
    T1DM(:) = Zero
    if (Spin) T1SDM(:) = Zero
    do q=1,NAC
      do p=1,NAC
        T1DM(p+(q-1)*NAC) = TS(q+(p-1)*NAC)
        if (Spin) T1SDM(p+(q-1)*NAC) = TSS(q+(p-1)*NAC)
      end do
    end do
    if (CI_Verify('SIGMA')) then
      call mma_allocate(BL,nVec,Label='BL')
      call mma_allocate(KL,nVec,Label='KL')
      call SG_Reord(iState,CI_St(iState)%STSym,1,nVec,Bra_Vec,BL)
      call SG_Reord(iState,CI_St(iState)%STSym,1,nVec,Ket_Vec,KL)
      call mma_allocate(DTMP,max(nT1DM,NAC*NAC),Label='DTMP')
      call mma_allocate(DSTMP,max(nT1DM,NAC*NAC),Label='DSTMP')
      call Lucia_Util('Densi',CI_Vector=KL(:),RVec=BL(:))
      ! (Lucia computes no spin density for singlets)
      Dev(1) = maxval(abs(DTMP(1:NAC*NAC)-T1DM(1:NAC*NAC)))
      Dev(2) = Zero
      if (Spin .and. (SGS(iState)%iSpin > 1)) Dev(2) = maxval(abs(DSTMP(1:NAC*NAC)-T1SDM(1:NAC*NAC)))
      if (.not. T1DM_SG_Verified) then
        write(u6,'(1X,A,2ES10.3)') 'Mk_T1DM (SGCI): SGUGA T1DM, T1SDM vs Lucia, max. abs. differences',Dev(1:2)
        T1DM_SG_Verified = .true.
      end if
      if (any(Dev(1:2) > ThrD)) then
        write(u6,'(1X,A,2ES10.3)') 'SGUGA error in the transition densities (SGCI), max. abs. differences',Dev(1:2)
        call Abend()
      end if
      call mma_deallocate(DSTMP)
      call mma_deallocate(DTMP)
      call mma_deallocate(KL)
      call mma_deallocate(BL)
    end if
    call mma_deallocate(TSS)
    call mma_deallocate(TS)
  end subroutine SGCI_T1DM

End Subroutine Mk_T1DM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> One-particle spin density SD(iTri(p,q)) = <Psi|a+_pa a_qa - a+_pb a_qb|Psi>
!> (levels, triangular packing as sg_spin_pdm) of the CI vector CI in the
!> split-graph order, from the ASTA triplet operator (SG_Tpq_Psi): the
!> segment values of the triplet operator, no two-particle intermediates.
subroutine SG_Spin_Density_ASTA(iState,LSM,n,CI,nTri,SD)
use sguga, only: SG_Tpq_Psi
use Index_Functions, only: iTri
integer(kind=iwp), intent(in) :: iState, LSM, n, nTri
real(kind=wp), intent(in) :: CI(n)
real(kind=wp), intent(out) :: SD(nTri)
integer(kind=iwp) :: p, q
real(kind=wp), allocatable :: Sgm(:)
call mma_allocate(Sgm,n,Label='SD_Sgm')
SD(:) = Zero
do p=1,SGS(iState)%nLev
  do q=1,p
    if (SGS(iState)%ISm(p) /= SGS(iState)%ISm(q)) cycle
    Sgm(:) = Zero
    call SG_Tpq_Psi(iState,p,q,One,LSM,CI,Sgm)
    SD(iTri(p,q)) = dot_product(CI,Sgm)
  end do
end do
call mma_deallocate(Sgm)
end subroutine SG_Spin_Density_ASTA

Subroutine Mk_pdms(CIVec,nCIVEC,D,SD,P,PA,nD,nP)
use lucia_data, only: DStmp, Dtmp, PAtmp, PTmp
use rasscf_global, only: NACPAR, NACPR2

 implicit none
 integer(kind=iwp), intent(in) :: nCIVEC
 real(kind=wp), intent(inout) :: CIVEC(nCIVEC)
 real(kind=wp), intent(out), optional :: D(:), SD(:), P(:), PA(:)
 integer(kind=iwp), intent(in) :: nD, nP

 real(kind=wp), allocatable :: D_loc(:), SD_loc(:), P_loc(:), PA_loc(:)
 real(kind=wp), allocatable :: D_FAROALD(:,:)
 real(kind=wp), allocatable :: SD_FAROALD(:,:)
 real(kind=wp), allocatable :: Faroald_Psi(:,:)
 real(kind=wp), allocatable :: P_Faroald(:,:,:,:)
 real(kind=wp), allocatable :: CIV(:)
 integer(kind=iwp), parameter :: iState=1

real(kind=wp) :: Check_D1, Check_P, Check_PA
real(kind=wp), allocatable :: D_Sguga(:)
real(kind=wp), allocatable :: P_Sguga(:), PA_sguga(:)

 ! SGCI: the densities from SGUGA in its basis (D, SD in the triangular
 ! packing of the levels, P and PA as sg_two_pdm); under SIGMA compared
 ! with those of Lucia for the vector reordered to the Lucia basis.
 if (SG_CI()) then
   call SGCI_Densities()
   return
 end if

 If (DoFaro) Then
   call mma_allocate(D_loc,nD,Label='D_loc')
   call mma_allocate(SD_loc,nD,Label='SD_loc')
   call mma_allocate(P_loc,nP,Label='P_loc')
   call mma_allocate(PA_loc,nP,Label='PA_loc')

   Call mma_allocate(D_Faroald,NAC,NAC)
   Call mma_allocate(SD_Faroald,NAC,NAC)
   Call mma_allocate(P_Faroald,NAC,NAC,NAC,NAC)
   If (ras_active) Then
     ! RAS: in the RAS string classes (faroald_ras)
     Call mma_allocate(Faroald_Psi,max(1,n_ras),1,Label='Psi')
     call CITRANS_CSF2RAS(CIVEC,Faroald_Psi(:,1))
     Call RAS_Trans_One_PDM(Faroald_Psi(:,1),Faroald_Psi(:,1),D_Faroald,SD_Faroald)
     Call RAS_Trans_Two_PDM(Faroald_Psi(:,1),Faroald_Psi(:,1),P_Faroald)
   Else
     Call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Psi')

     ! (the CI vector is in the configuration order of Faroald, citrans)
     Faroald_Psi(:,:)=Zero
     call CITRANS_CSF2SD(CIVEC,Faroald_PSI)

     Call One_pdm(Faroald_Psi,D_Faroald,SD_Faroald)
     Call two_pdm(Faroald_psi,P_Faroald)
   End If
   Call Fold2(1,[NAC],D_faroald,D_loc)
   Call Fold2(1,[NAC],SD_faroald,SD_loc)
   Call Fold_Two_pdm(P_Faroald,P_loc,PA_loc)

   Call mma_deallocate(Faroald_Psi)
   Call mma_deallocate(P_faroald)
   Call mma_deallocate(D_faroald)
   Call mma_deallocate(SD_faroald)

   If (Present(D)) D(1:nD)=D_loc(1:nD)
   If (Present(SD)) SD(1:nD)=SD_loc(1:nD)
   If (Present(P)) P(1:nP)=P_loc(1:nP)
   If (Present(PA)) PA(1:nP)=PA_loc(1:nP)

   If (CI_Verify('SIGMA')) call Faroald_Check_Densities()

   call mma_deallocate(D_loc)
   call mma_deallocate(SD_loc)
   call mma_deallocate(P_loc)
   call mma_deallocate(PA_loc)
 Else
   Call mma_allocate(DTMP,nD,Label='DTMP')
   Call mma_allocate(DSTMP,nD,Label='DSTMP')
   Call mma_allocate(PTMP,nP,Label='PTMP')
   Call mma_allocate(PATMP,nP,Label='PATMP')
   call Lucia_Util('Densi',CI_Vector=CIVEC)
   If (Present(D)) D(1:nD)=DTMP(1:nD)
   If (Present(SD)) SD(1:nD)=DSTMP(1:nD)
   If (Present(P)) P(1:nP)=PTMP(1:nP)
   If (Present(PA)) PA(1:nP)=PATMP(1:nP)
   Call mma_deallocate(DTMP)
   Call mma_deallocate(DSTMP)
   Call mma_deallocate(PTMP)
   Call mma_deallocate(PATMP)
 End If

! CI development (MOLCAS_CI_VERIFY=SIGMA): verify the SGUGA density matrices
! against the CI back end. D, P and PA are optional, so only those present
! are compared.
        If (CI_Verify('SIGMA') .and. (.NOT.iDoGAS) .and. (.NOT.DoFaro)) Then

          Call mma_allocate(CIV,nCIVEC,Label='CIV')
          call CI_Basis_To_File(iState,CI_St(iState)%STSym,nCIVEC,CIVEC,CIV)

!         Test the one-particle density matrix
          If (Present(D)) Then
          Check_D1=CheckSum(D,NACPAR)
          Call mma_allocate(D_sguga,NAC*(NAC+1)/2)

          call sg_one_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),CI_St(iState)%STSym,D_sguga,Size(D_sguga))
          If (ABS(CheckSum(D_sguga,NACPAR)-Check_D1)/SIZE(D_sguga)>1.0e-12_wp) Then
             Check_D1=CheckSum(D_sguga,NACPAR)
             Write (u6,*) 'SGUGA error in D1Mat'
             Call Abend()
          End If
          Call mma_deallocate(D_sguga)
          End If

!         Test the one-particle spin-density matrix
!         This option is not yet developed for the SGUGA code. To come...

!         Test the symmetric two-particle density matrix.
          If (Present(P) .and. Present(PA)) Then
          Check_P=CheckSum(P,NACPR2)
          Check_PA=CheckSum(PA,NACPR2)

          Call mma_allocate(P_sguga,NACPR2,Label='P')
          Call mma_allocate(PA_sguga,NACPR2,Label='PA')

          call SG_Two_PDM_G(iState,CIV,SIZE(CIV),CI_St(iState)%STSym,P_sguga,PA_sguga,NACPAR*(NACPAR+1)/2)

          If (ABS(CheckSum(P_sguga,NACPR2)-Check_P)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_P=CheckSum(P_sguga,NACPR2)
             Write (u6,*) 'SGUGA error in P'
             Call Abend()
          End If

          If (ABS(CheckSum(PA_sguga,NACPR2)-Check_PA)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_PA=CheckSum(PA_sguga,NACPR2)
             Write (6,*) 'SGUGA error in PA'
             Call Abend()
          End If

          Call mma_deallocate(PA_sguga)
          Call mma_deallocate(P_sguga)
          End If
          Call mma_deallocate(CIV)

        END IF

contains

  ! FARO (CI development, SIGMA): the Faroald densities D, SD, P, PA against
  ! Lucia (vector in the Lucia basis; no spin density for singlets) and
  ! SGUGA (vector in the split-graph order), element by element.
  ! CI development (SIGMA): the spin density from the ASTA triplet operator
  ! (SG_Spin_Density_ASTA, vector CIS in the split-graph order) against SDRef.
  subroutine Check_SD_ASTA(CIS_,SDRef,Text)
    real(kind=wp), intent(in) :: CIS_(:), SDRef(:)
    character(len=*), intent(in) :: Text
    real(kind=wp), allocatable :: SDA(:)
    logical(kind=iwp) :: Done
    call mma_allocate(SDA,NACPAR,Label='SDA')
    call SG_Spin_Density_ASTA(iState,CI_St(iState)%STSym,nCIVEC,CIS_,NACPAR,SDA)
    Done = SD_ASTA_Verified
    call Report_Dev(Text,maxval(abs(SDA(:)-SDRef(1:NACPAR))),1.0e-10_wp,Done)
    SD_ASTA_Verified = .true.
    call mma_deallocate(SDA)
  end subroutine Check_SD_ASTA

  subroutine Faroald_Check_Densities()
    real(kind=wp), allocatable :: CIL(:), DS(:), SDS(:), PS(:), PAS(:)
    real(kind=wp) :: Dev(4)
    logical(kind=iwp) :: Done
    real(kind=wp), parameter :: ThrD = 1.0e-10_wp

    call mma_allocate(CIL,nCIVEC,Label='CIL')
    call CI_Basis_To_Lucia(iState,CI_St(iState)%STSym,nCIVEC,CIVEC,CIL)
    call mma_allocate(DTMP,nD,Label='DTMP')
    call mma_allocate(DSTMP,nD,Label='DSTMP')
    call mma_allocate(PTMP,nP,Label='PTMP')
    call mma_allocate(PATMP,nP,Label='PATMP')
    call Lucia_Util('Densi',CI_Vector=CIL)
    Dev(1) = maxval(abs(DTMP(1:NACPAR)-D_loc(1:NACPAR)))
    Dev(2) = Zero
    if (CI_St(iState)%iSpin > 1) Dev(2) = maxval(abs(DSTMP(1:NACPAR)-SD_loc(1:NACPAR)))
    Dev(3) = maxval(abs(PTMP(1:nP)-P_loc(1:nP)))
    Dev(4) = maxval(abs(PATMP(1:nP)-PA_loc(1:nP)))
    Done = PDM_Faro_Verified
    call Report_Dev('Mk_pdms (FARO): Faroald D, SD, P, PA vs Lucia',maxval(Dev),ThrD,Done)
    call mma_deallocate(PATMP)
    call mma_deallocate(PTMP)
    call mma_deallocate(DSTMP)
    call mma_deallocate(DTMP)
    if (SG_Avail()) then
      call CI_Basis_To_File(iState,CI_St(iState)%STSym,nCIVEC,CIVEC,CIL)
      call mma_allocate(DS,NACPAR,Label='DS')
      call mma_allocate(SDS,NACPAR,Label='SDS')
      call mma_allocate(PS,NACPR2,Label='PS')
      call mma_allocate(PAS,NACPR2,Label='PAS')
      call sg_one_pdm(SGS(iState),CIS(iState),EXS(iState),CIL,nCIVEC,CI_St(iState)%STSym,DS,NACPAR)
      call SG_Spin_PDM_G(iState,CIL,nCIVEC,CI_St(iState)%STSym,DS,SDS,NACPAR)
      call SG_Two_PDM_G(iState,CIL,nCIVEC,CI_St(iState)%STSym,PS,PAS,NACPR2)
      Dev(1) = maxval(abs(DS(:)-D_loc(1:NACPAR)))
      Dev(2) = maxval(abs(SDS(:)-SD_loc(1:NACPAR)))
      call Check_SD_ASTA(CIL,SD_loc(1:NACPAR),'Mk_pdms (FARO): spin density, ASTA T_pq vs Faroald')
      Dev(3) = maxval(abs(PS(1:nP)-P_loc(1:nP)))
      Dev(4) = maxval(abs(PAS(1:nP)-PA_loc(1:nP)))
      Done = PDM_Faro_Verified
      call Report_Dev('Mk_pdms (FARO): Faroald D, SD, P, PA vs SGUGA',maxval(Dev),ThrD,Done)
      call mma_deallocate(PAS)
      call mma_deallocate(PS)
      call mma_deallocate(SDS)
      call mma_deallocate(DS)
    end if
    PDM_Faro_Verified = .true.
    call mma_deallocate(CIL)
  end subroutine Faroald_Check_Densities

  subroutine SGCI_Densities()
    real(kind=wp), allocatable :: DS(:), SDS(:), PS(:), PAS(:), CIL(:)
    real(kind=wp) :: Dev(4)
    real(kind=wp), parameter :: ThrD = 1.0e-10_wp

    call mma_allocate(DS,NACPAR,Label='DS')
    call mma_allocate(SDS,NACPAR,Label='SDS')
    call mma_allocate(PS,NACPR2,Label='PS')
    call mma_allocate(PAS,NACPR2,Label='PAS')
    call sg_one_pdm(SGS(iState),CIS(iState),EXS(iState),CIVEC,nCIVEC,CI_St(iState)%STSym,DS,NACPAR)
    SDS(:) = Zero
    if (Present(SD)) call SG_Spin_PDM_G(iState,CIVEC,nCIVEC,CI_St(iState)%STSym,DS,SDS,NACPAR)
    PS(:) = Zero
    PAS(:) = Zero
    if (Present(P) .or. Present(PA)) call SG_Two_PDM_G(iState,CIVEC,nCIVEC,CI_St(iState)%STSym,PS,PAS,NACPR2)
    if (Present(D)) then
      D(1:nD) = Zero
      D(1:NACPAR) = DS(:)
    end if
    if (Present(SD)) then
      SD(1:nD) = Zero
      SD(1:NACPAR) = SDS(:)
    end if
    if (Present(P)) P(1:nP) = PS(1:nP)
    if (Present(PA)) PA(1:nP) = PAS(1:nP)

    if (CI_Verify('SIGMA')) then
      ! Lucia densities of the vector in the Lucia basis
      call mma_allocate(CIL,nCIVEC,Label='CIL')
      call SG_Reord(iState,CI_St(iState)%STSym,1,nCIVEC,CIVEC,CIL)
      call mma_allocate(DTMP,nD,Label='DTMP')
      call mma_allocate(DSTMP,nD,Label='DSTMP')
      call mma_allocate(PTMP,nP,Label='PTMP')
      call mma_allocate(PATMP,nP,Label='PATMP')
      call Lucia_Util('Densi',CI_Vector=CIL)
      ! Only the triangles are defined (DTMP and DSTMP are not initialized
      ! behind them), and Lucia computes no spin density for singlets.
      Dev(1) = maxval(abs(DTMP(1:NACPAR)-DS(:)))
      Dev(2) = Zero
      if (Present(SD) .and. (SGS(iState)%iSpin > 1)) Dev(2) = maxval(abs(DSTMP(1:NACPAR)-SDS(:)))
      if (Present(SD)) call Check_SD_ASTA(CIVEC,SDS,'Mk_pdms (SGCI): spin density, ASTA T_pq vs sg_spin_pdm')
      Dev(3) = Zero
      Dev(4) = Zero
      if (Present(P) .or. Present(PA)) then
        Dev(3) = maxval(abs(PTMP(1:nP)-PS(1:nP)))
        Dev(4) = maxval(abs(PATMP(1:nP)-PAS(1:nP)))
      end if
      if (.not. PDM_SG_Verified) then
        write(u6,'(1X,A,4ES10.3)') 'Mk_pdms (SGCI): SGUGA D, SD, P, PA vs Lucia, max. abs. differences',Dev(:)
        PDM_SG_Verified = .true.
      end if
      if (any(Dev(:) > ThrD)) then
        write(u6,'(1X,A,4ES10.3)') 'SGUGA error in the density matrices (SGCI), max. abs. differences',Dev(:)
        call Abend()
      end if
      call mma_deallocate(PATMP)
      call mma_deallocate(PTMP)
      call mma_deallocate(DSTMP)
      call mma_deallocate(DTMP)
      call mma_deallocate(CIL)
    end if

    call mma_deallocate(PAS)
    call mma_deallocate(PS)
    call mma_deallocate(SDS)
    call mma_deallocate(DS)
  end subroutine SGCI_Densities

 End Subroutine Mk_pdms

!***********************************************************************************************************************************
!***********************************************************************************************************************************

 Function Checksum(A,nA)
 real(kind=wp) :: Checksum
 integer(kind=iwp), intent(in):: nA
 real(kind=wp), intent(in):: A(nA)
 integer(kind=iwp) :: i
 Checksum=0.0_wp
 Do i = 1, nA
    Checksum = Checksum + Abs(A(i))/real(i,kind=wp)
 End Do
 End Function Checksum

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine CI_Timer(Opt)
use lucia_data, only: TDENSI, TSIGMA

character(Len=*) Opt

Select case (Opt(1:3))
Case ('Ini')
   !SVC: lucia timers
   tsigma(:) = Zero
   tdensi(:) = Zero
Case ('Pri')
  If (DoFaro) Then
  Else
    call PrtTim()
  End If
Case Default
  Write (u6,*) 'CI_Timer called with illegal option:',Opt
  Call Abend()
End Select

end subroutine CI_Timer

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> The CI space of RASSCF (Space_CI) from the GAS description of the input
!> (NGAS, NGSSH, IGSOCCX; CAS and RAS are converted to it by proc_inp) and
!> the HEXS/DEXS lists. With MOLCAS_CI_VERIFY (any level) printed and
!> checked against the RAS description (nRsPrt, nRas, nRasEl, NSM) and the
!> HEXS/DEXS thresholds used by the CI interface.
subroutine CI_Space_Setup()
use general_data, only: IGSOCCX, ISPIN, NACTEL, NGAS, NGSSH, nElec3, nHole1, nRs1, nRs2, nRs3, nSym
use CI_Info, only: Space_Inp
use CI_States, only: RAS_Partitions
use Molcas, only: MxSym
integer(kind=iwp) :: g, k, i, iLev, iSym, iPart(MxGAS), nBelow, nAbove, MinR, MaxR, MinEO(MxGAS), nRsPrt, nRas(MxSym,3), &
                     nRasEl(3)
integer(kind=iwp), allocatable :: NMx(:)
logical(kind=iwp) :: OK

call Space_From_GAS(Space_CI,nSym,NACTEL,ISPIN,NGAS,NGSSH(1:NGAS,1:nSym),IGSOCCX(1:NGAS,1:2))
! partition of each input space (0: no orbitals, dropped)
k = 0
do g=1,NGAS
  iPart(g) = 0
  if (sum(NGSSH(g,1:nSym)) == 0) cycle
  k = k+1
  iPart(g) = k
end do
! HEXS (maximum occupation) and DEXS (maximum - 1) removals
if (Space_Removals(Space_Inp)) then
  call mma_allocate(NMx,max(1,Space_CI%nPart),Label='NMx')
  call Space_Max_Occ(Space_CI,NMx)
  do i=1,Space_Inp%nHEXS
    g = Space_Inp%iHEXS(i)
    if (iPart(g) > 0) call Space_Remove(Space_CI,iPart(g),NMx(iPart(g)))
  end do
  do i=1,Space_Inp%nDEXS
    g = Space_Inp%iDEXS(i)
    if (iPart(g) > 0) call Space_Remove(Space_CI,iPart(g),NMx(iPart(g))-1)
  end do
  call mma_deallocate(NMx)
end if
if (.not. CI_Verify('ANY')) return

call Space_Print(Space_CI)
! GAS: no RAS description to compare with
if (iDoGAS) return
! the RAS description (as CI_Setup_RAS builds it), with its empty
! partitions dropped (their lower
! bound merged into the previous one): the same partitions, the same lower
! bounds (after tightening to the physical ones), no upper bounds, the same
! levels
call RAS_Partitions(nSym,NACTEL,nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,.false.,nRsPrt,nRas,nRasEl)
OK = (Space_CI%nLev == NAC)
k = 0
nBelow = 0
do i=1,nRsPrt
  if (sum(nRas(1:nSym,i)) == 0) then
    if (k > 0) MinEO(k) = max(MinEO(k),nRasEl(i))
    cycle
  end if
  k = k+1
  if (k > Space_CI%nPart) then
    OK = .false.
    exit
  end if
  OK = OK .and. all(Space_CI%nOrb(1:nSym,k) == nRas(1:nSym,i))
  MinEO(k) = nRasEl(i)
end do
OK = OK .and. (k == Space_CI%nPart)
if (OK) then
  do k=1,Space_CI%nPart
    nBelow = nBelow+sum(Space_CI%nOrb(1:nSym,k))
    nAbove = sum(Space_CI%nOrb(1:nSym,k+1:Space_CI%nPart))
    MinR = max(MinEO(k),0,NACTEL-2*nAbove)
    MaxR = min(NACTEL,2*nBelow)
    if (k == Space_CI%nPart) then
      MinR = NACTEL
      MaxR = NACTEL
    end if
    OK = OK .and. (Space_CI%MinEl(k) == MinR) .and. (Space_CI%MaxEl(k) == MaxR)
  end do
  ! the levels: partition by partition, within each by symmetry
  iLev = 0
  do i=1,nRsPrt
    do iSym=1,nSym
      OK = OK .and. all(Space_CI%LevSym(iLev+1:iLev+nRas(iSym,i)) == iSym)
      iLev = iLev+nRas(iSym,i)
    end do
  end do
end if
if (.not. OK) then
  write(u6,'(1X,A)') 'CI_Space_Setup: the CI space differs from the RAS description (nRsPrt, nRas, nRasEl)'
  call Abend()
end if
write(u6,'(1X,A)') 'CI_Space_Setup: CI space verified against the RAS description'
end subroutine CI_Space_Setup

subroutine CI_Initialize(DO_SGUGA,DO_FAROALD)
! (the input values: the CI state is set up here)
use general_data, only: ISPIN, NACTEL, STSYM
use sguga, only: SG_Test_EpqErs_Ref
use faroald, only: nela, nelb
use second_quantization, only: binom_coef
use general_data, only: STSym, nSym
use Molcas, only: MxSym

implicit none
logical(kind=iwp), intent(in) :: DO_SGUGA,DO_FAROALD
integer(kind=iwp), parameter :: iState=1
integer(kind=iwp) :: nAct, k, iSym, iLev, nPrt, nCSF_Lucia
integer(kind=iwp), allocatable :: nOrbP(:)
logical(kind=iwp) :: Full
integer(kind=iwp), allocatable :: OrbSym(:), OrbRas(:), nOff(:)

nCSF_Lucia = 0
! SGCI does not need Lucia; it is set up only as the reference of
! MOLCAS_CI_VERIFY=SIGMA.
Lucia_Active = .not. (((Do_SGUGA .and. SG_CI()) .or. Do_Faroald) .and. (.not. CI_Verify('SIGMA')))

! the CI space in the common description of the back ends
call CI_Space_Setup()

if (Lucia_Active) then
  ! Construct the determinant tables
  ! Initialize LUCIA and determinant control
  call Lucia_Util('Ini')

  ! to get number of CSFs for GAS
  ! and number of determinants to store (Lucia's dimensions, CI_St set by
  ! lucia2molcas; Lucia's CSF count of the CI symmetry kept for the check
  ! against the split graph, which redefines CI_St%nCSF below)
  nconf = MaxVal(CI_St(iState)%nCSF(1:mxsym))
  nDet = MaxVal(CI_St(iState)%nDet(1:mxsym))
  nCSF_Lucia = CI_St(iState)%nCSF(STSym)
end if

! Initiate the SGUGA environment conditional to all flags
If (Do_SGUGA) Then
   if (CI_Verify('SIGMA')) then
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
      write(u6,'(1X,A)') 'CI back end cross-checked against SGUGA (MOLCAS_CI_VERIFY=SIGMA).'
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
   end if
   call SG_Setup_RASSCF()
   ! GAS: the split graph against Lucia's CSF count of the CI symmetry (the
   ! CI is Lucia's; with GAS Lucia counts the CI symmetry only)
   if (iDoGAS .and. Lucia_Active) then
     if (CIS(iState)%nCSF(STSym) /= nCSF_Lucia) then
       write(u6,'(1X,A,2I10)') 'CI_Initialize: GAS split graph and Lucia, CSFs of the CI symmetry',CIS(iState)%nCSF(STSym), &
                               nCSF_Lucia
       call Abend()
     end if
     if (CI_Verify('ANY')) write(u6,'(1X,A,I0,A)') 'CI_Initialize: GAS split graph verified against Lucia, ', &
                                                   CIS(iState)%nCSF(STSym),' CSFs of the CI symmetry'
   end if
   if (.not. Lucia_Active) then
     ! the CSF counts from SGUGA; no determinants (the Davidson scratch for
     ! the determinant expansion of the other back ends is not needed)
     ! (CI_St%nCSF from the split graph, CI_Setup_SG)
     CI_St(iState)%nCSF(:) = 0
     CI_St(iState)%nCSF(1:nSym) = CIS(iState)%nCSF(1:nSym)
     CI_St(iState)%nDet(:) = CI_St(iState)%nCSF(:)
     nDet = nConf
     write(u6,'(1X,A)') 'SGCI/FARO: Lucia is not initialized.'
   end if
   ! Exhaustive verification of the SGUGA coupling-coefficient codes,
   ! only on request (used by the test/csf_sd groups), see CI_Verify. With
   ! GAS only the single operators (E_pq, T_pq), which are exact in any
   ! restricted graph: the products of the self-test sum over the
   ! intermediate CSFs of the graph itself, which is exact for CAS/RAS only
   ! (the CI interface evaluates them for GAS in the relaxed space, checked
   ! under SIGMA against Lucia and Faroald).
   if (CI_Verify('YES') .and. (.not. iDoGAS)) then
      call SG_Test_EpqErs_Ref(iState,STSym,sg_h_psi,sg_h_psi_ras)
   else if (CI_Verify('YES') .or. CI_Verify('EPQ')) then
      call SG_Test_EpqErs_Ref(iState,STSym,sg_h_psi,sg_h_psi_ras,Epq_Only=.true.)
   end if
End If

! faroald initializations
if (DO_FAROALD) then
  write(u6,'(1X,A)') '**EXPERIMENTAL**'
  write(u6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
  write(u6,'(1X,A)') '**EXPERIMENTAL**'
  ! the orbitals of Faroald = the levels of the split graph: by partition
  ! (RAS/GAS), within each by symmetry (Space_CI%LevSym);
  ! Faro_L2O maps them to the symmetry-blocked order of RASSCF (within one
  ! symmetry the orbitals keep the order of their levels)
  nAct = sum(NASH(1:nSym))
  call mma_allocate(OrbSym,max(1,nAct),Label='OrbSym')
  call mma_allocate(OrbRas,max(1,nAct),Label='OrbRas')
  call mma_deallocate(Faro_L2O,safe='*')
  call mma_allocate(Faro_L2O,max(1,nAct),Label='Faro_L2O')
  iLev = Space_CI%nLev
  if (iLev == nAct) then
    OrbSym(1:nAct) = Space_CI%LevSym(1:nAct)
    OrbRas(1:nAct) = Space_CI%LevPart(1:nAct)
  end if
  nPrt = Space_CI%nPart
  if ((iLev /= nAct) .or. (NAC /= nAct)) then
    write(u6,*) 'CI_Initialize: inconsistent active levels for Faroald',iLev,nAct,NAC
    call Abend()
  end if
  call mma_allocate(nOff,nSym,Label='nOff')
  nOff(1) = 0
  do iSym=2,nSym
    nOff(iSym) = nOff(iSym-1)+NASH(iSym-1)
  end do
  do iLev=1,nAct
    nOff(OrbSym(iLev)) = nOff(OrbSym(iLev))+1
    Faro_L2O(iLev) = nOff(OrbSym(iLev))
  end do
  call mma_deallocate(nOff)
  ! the tables of the full determinant space: for CAS; with RAS/GAS only
  ! as the reference of MOLCAS_CI_VERIFY=SIGMA, for small CAS spaces
  Full = (nPrt == 1)
  if ((.not. Full) .and. CI_Verify('SIGMA')) &
    Full = real(binom_coef((NACTEL+ISPIN-1)/2,nAct),kind=wp)*real(binom_coef((NACTEL-ISPIN+1)/2,nAct),kind=wp) <= 1.0e6_wp
  call FAROALD_INIT(NACTEL,nAct,ISPIN,Full)
  call CITRANS_INIT(NACTEL,nAct,ISPIN,nSym,OrbSym,STSym,nPrt,OrbRas,Space_CI%MinEl(1:nPrt),Space_CI%MaxEl(1:nPrt))
  ! RAS/GAS: the sigma vector in the string classes (faroald_ras)
  if (nPrt > 1) then
    call mma_allocate(nOrbP,nPrt,Label='nOrbP')
    do k=1,nPrt
      nOrbP(k) = count(OrbRas(1:nAct) == k)
    end do
    call faroald_ras_init(nela,nelb,nAct,nPrt,nOrbP,Space_CI%MinEl(1:nPrt),Space_CI%MaxEl(1:nPrt),nSym,OrbSym(1:nAct),STSym, &
                          ISPIN == 1)
    call mma_deallocate(nOrbP)
    write(u6,'(1X,A,I0,A,I0,A)') 'Faroald: string-driven RAS/GAS, ',nPrt,' partitions, ',n_ras, &
                                 ' determinants of the CI symmetry'
  end if
  call mma_deallocate(OrbSym)
  call mma_deallocate(OrbRas)
  if ((citrans_ncsf() /= nConf) .and. (nActEl > 0)) then
    write(u6,*) 'CI_Initialize: number of CSFs, Faroald and RASSCF:',citrans_ncsf(),nConf
    call Abend()
  end if
end if

! HEXS/DEXS: the kept CSFs
call HEXS_Setup(iState,Do_SGUGA)

End subroutine CI_Initialize

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine CI_Close(iState)
integer(kind=iwp), Intent(in):: iState

! All of these are safe to call even of their corresponding environment have not been activated.

if (Lucia_Active) Call Lucia_Util('Close')

call faroald_free()
call faroald_ras_free()
call citrans_free()
call mma_deallocate(Faro_L2O,safe='*')

if (CI_St(iState)%iRel > 0) then
  call SG_Free(CI_St(iState)%iRel)
  CI_St(iState)%iRel = 0
end if
call SG_Rel_Free()
Call SG_free(iState)
call SG_Free_Cnf_List(Cnf_SG)
call mma_deallocate(HEXS_Keep,safe='*')
call mma_deallocate(HEXS_Keep_SG,safe='*')
call Space_Free(Space_CI)

end subroutine CI_Close

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Diagonal of the CI Hamiltonian in the CSF basis as used by the CI solver
!> (preconditioner, selection of the explicit subspace in PHPCSF): for each
!> configuration the average of its determinant diagonals, assigned to all
!> of its CSFs (Lucia, CSDIAG). With DoFaro the same average from Faroald
!> (hdiag), compared with Lucia and then used.
!>
!> With MOLCAS_CI_VERIFY=SIGMA (single RAS partition, no GAS) the exact
!> diagonal <I|H|I> is computed by SGUGA (sg_hdiag_fast, for small spaces
!> checked against sg_hdiag_ref) and compared (Verify_IHI): exactly for the
!> single-determinant configurations, where the average is exact, and
!> otherwise by the configurations that the selection of PHPCSF would pick
!> with either diagonal. The latter is reported, not enforced, since the
!> average is an approximation.
subroutine Mk_IHI(CSFDIA,nConf,TU,nTU,TUVX,nTUVX)
use Lucia_Interface, only: Lucia_Util
use Lucia_data, only: SDREO
use lucia_data, only: MINOP, NCNFTP, NCSFTP, NDTFTP, NTYP
use Constants, only: One

integer(kind=iwp), intent(in) :: nConf, nTU, nTUVX
real(kind=wp), intent(out):: CSFDIA(nConf)
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

real(kind=wp), allocatable :: DDIA(:)
integer(kind=iwp):: IPRINT=0

real(kind=wp), parameter :: ThrDia = 1.0e-10_wp
! the CI state of RASSCF (SGS, CIS, EXS)
integer(kind=iwp), parameter :: iState = 1

! SGCI: the exact diagonal in the SGUGA basis, no Lucia
if (SG_CI()) then
  call SGCI_Diag()
  return
end if
! Faroald (FARO): the configuration average in its own basis
if (DoFaro) then
  call Faroald_Diag()
  return
end if

! COMPUTE CI DIAGONAL IN DETERMINANT BASIS

call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)

call mma_allocate(DDIA,NDET,label='DETDIA')

call get_diag(DDIA,ndet)

! Todo: ndet comes from general_data!! Check this!
! TRANSFORM CI DIAGONAL FROM DET TO CSF BASIS
call CSDIAG(NCONF,ndet,CSFDIA,DDIA,NCNFTP(1,CI_St(iState)%STSym),NTYP,SDREO,NDTFTP,NCSFTP,IPRINT)

! DEALLOCATE LOCAL MEMORY
call mma_deallocate(DDIA)

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the exact diagonal.
if (CI_Verify('SIGMA') .and. SG_Avail()) call Verify_IHI()

! CI development (MOLCAS_CI_VERIFY=HDIAG): use the exact diagonal instead.
if (CI_Verify('HDIAG') .and. SG_Avail()) then
  call Exact_Diag(CSFDIA)
  write(u6,'(1X,A)') 'Mk_IHI: exact SGUGA diagonal used (MOLCAS_CI_VERIFY=HDIAG)'
end if

contains

  ! Faroald (FARO): the CI-solver diagonal in its own basis (citrans), the
  ! configuration average of Lucia: for each configuration the average of
  ! the diagonal elements of its determinants (M_S = S), Faroald_HDet. Under
  ! SIGMA compared with Lucia (the diagonal reordered to the Lucia basis).
  subroutine Faroald_Diag()
    use citrans, only: citrans_conf_dets
    integer(kind=iwp) :: k, k1, nc, nd, i
    integer(kind=iwp), allocatable :: da(:), db(:)
    real(kind=wp) :: E, MaxDev
    real(kind=wp), allocatable :: cf(:,:), DD(:), DL(:), Idx(:), Map(:)

    call Faroald_Ints(TU,TUVX)
    k = 1
    do while (k <= nConf)
      call citrans_conf_dets(k,k1,nc,nd,da,db,cf)
      E = Zero
      do i=1,nd
        E = E+Faroald_HDet(da(i),db(i),da(i),db(i))
      end do
      CSFDIA(k1:k1+nc-1) = E/real(max(nd,1),kind=wp)
      k = k1+nc
      call mma_deallocate(da)
      call mma_deallocate(db)
      call mma_deallocate(cf)
    end do
    call HEXS_Proj(HEXS_Keep,CSFDIA)
    if (CI_Verify('SIGMA')) then
      call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
      call mma_allocate(DD,NDET,label='DETDIA')
      call get_diag(DD,ndet)
      call mma_allocate(DL,nConf,Label='DL')
      call CSDIAG(NCONF,ndet,DL,DD,NCNFTP(1,CI_St(iState)%STSym),NTYP,SDREO,NDTFTP,NCSFTP,IPRINT)
      ! the permutation Faroald -> Lucia basis
      call mma_allocate(Idx,nConf,Label='Idx')
      call mma_allocate(Map,nConf,Label='Map')
      do i=1,nConf
        Idx(i) = real(i,kind=wp)
      end do
      call CI_Basis_To_Lucia(iState,CI_St(iState)%STSym,nConf,Idx,Map)
      MaxDev = Zero
      do i=1,nConf
        MaxDev = max(MaxDev,abs(CSFDIA(nint(abs(Map(i))))-DL(i)))
      end do
      if (.not. IHI_Avg_Verified) &
        write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI (Faroald): configuration-average diagonal verified against Lucia, ',nConf, &
                                       ' CSFs, max. error ',MaxDev
      IHI_Avg_Verified = .true.
      if (MaxDev > ThrDia*max(One,maxval(abs(DL(:))))) then
        write(u6,*) 'Mk_IHI: the Faroald configuration-average diagonal differs from Lucia'
        call Abend()
      end if
      call mma_deallocate(Map)
      call mma_deallocate(Idx)
      call mma_deallocate(DL)
      call mma_deallocate(DD)
      if (SG_Avail()) call Faroald_Diag_SG()
    end if
  end subroutine Faroald_Diag

  ! The Faroald configuration-average diagonal against that of SGUGA
  ! (SG_Cnf_Avg_Diag, split-graph order of the files).
  subroutine Faroald_Diag_SG()
    use sguga, only: SG_Mk_Cnf_List, SG_Cnf_Avg_Diag
    type(SGCnfList) :: Cnf
    integer(kind=iwp) :: i
    real(kind=wp) :: MaxDev
    real(kind=wp), allocatable :: DSG(:), Idx(:), Map(:)
    call SG_Mk_Cnf_List(iState,CI_St(iState)%STSym,Cnf)
    call mma_allocate(DSG,nConf,Label='DSG')
    call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,DSG)
    call HEXS_Proj(HEXS_Keep_SG,DSG)
    call mma_allocate(Idx,nConf,Label='Idx')
    call mma_allocate(Map,nConf,Label='Map')
    do i=1,nConf
      Idx(i) = real(i,kind=wp)
    end do
    call CI_Basis_To_File(iState,CI_St(iState)%STSym,nConf,Idx,Map)
    MaxDev = Zero
    do i=1,nConf
      MaxDev = max(MaxDev,abs(CSFDIA(nint(abs(Map(i))))-DSG(i)))
    end do
    call Report_Dev('Mk_IHI (FARO): Faroald configuration-average diagonal vs SGUGA',MaxDev, &
                    ThrDia*max(One,maxval(abs(DSG(:)))),IHI_Faro_SG_Verified)
    call mma_deallocate(Map)
    call mma_deallocate(Idx)
    call mma_deallocate(DSG)
    call SG_Free_Cnf_List(Cnf)
  end subroutine Faroald_Diag_SG

  ! SGCI: the CI-solver diagonal in the split-graph CSF basis: the
  ! configuration average of Lucia (SG_Cnf_Avg_Diag), so that SGCI follows
  ! the iterations of Lucia, or with MOLCAS_CI_VERIFY=HDIAG the exact
  ! <I|H|I> (sg_hdiag_fast). Under SIGMA the exact diagonal is checked
  ! against sg_hdiag_ref (small spaces) and the average against Lucia.
  subroutine SGCI_Diag()
    use sguga, only: SG_Cnf_Avg_Diag
    integer(kind=iwp), parameter :: MxRef = 500
    real(kind=wp) :: MaxDev
    real(kind=wp), allocatable :: DRef(:), DAvg(:), DAvgL(:), DL(:), DD(:)
    call SG_Cnf_Kept(CI_St(iState)%STSym)
    call mma_allocate(DAvg,nConf,Label='DAvg')
    call SG_Cnf_Avg_Diag(iState,Cnf_SG,TU,nTU,TUVX,nTUVX,DAvg)
    call HEXS_Proj(HEXS_Keep_SG,DAvg)
    if (CI_Verify('SIGMA')) then
      ! Lucia's configuration-average diagonal, in the Lucia basis
      call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
      call mma_allocate(DD,NDET,label='DETDIA')
      call get_diag(DD,ndet)
      call mma_allocate(DL,nConf,Label='DL')
      call CSDIAG(NCONF,ndet,DL,DD,NCNFTP(1,CI_St(iState)%STSym),NTYP,SDREO,NDTFTP,NCSFTP,IPRINT)
      call mma_allocate(DAvgL,nConf,Label='DAvgL')
      call Diag_to_Lucia(DAvg,DAvgL)
      MaxDev = maxval(abs(DAvgL(:)-DL(:)))
      if (.not. IHI_Avg_Verified) &
        write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI (SGCI): SGUGA configuration-average diagonal verified against Lucia, ',nConf, &
                                       ' CSFs, max. error ',MaxDev
      IHI_Avg_Verified = .true.
      if (MaxDev > ThrDia*max(One,maxval(abs(DL(:))))) then
        write(u6,*) 'Mk_IHI: the SGUGA configuration-average diagonal differs from Lucia'
        call Abend()
      end if
      call mma_deallocate(DAvgL)
      call mma_deallocate(DL)
      call mma_deallocate(DD)
    end if
    if (.not. CI_Verify('HDIAG')) then
      CSFDIA(:) = DAvg(:)
      call mma_deallocate(DAvg)
      if (.not. CI_Verify('SIGMA')) return
    else
      call mma_deallocate(DAvg)
    end if
    ! the exact diagonal (HDIAG, and its check under SIGMA)
    call mma_allocate(DAvg,nConf,Label='DExact')
    call SG_HDiag_G('FAST',iState,nConf,CI_St(iState)%STSym,TUVX,nTUVX,TU,nTU,DAvg)
    if (CI_Verify('HDIAG')) then
      CSFDIA(:) = DAvg(:)
      call HEXS_Proj(HEXS_Keep_SG,CSFDIA)
      write(u6,'(1X,A)') 'Mk_IHI (SGCI): exact SGUGA diagonal used (MOLCAS_CI_VERIFY=HDIAG)'
    end if
    if (CI_Verify('SIGMA') .and. (nConf <= MxRef)) then
      call mma_allocate(DRef,nConf,Label='DRef')
      call SG_HDiag_G('REF',iState,nConf,CI_St(iState)%STSym,TUVX,nTUVX,TU,nTU,DRef)
      MaxDev = maxval(abs(DAvg(:)-DRef(:)))
      if (.not. IHI_SG_Verified) &
        write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI (SGCI): sg_hdiag_fast verified against sg_hdiag_ref, ',nConf,' CSFs, max. error ', &
                                       MaxDev
      IHI_SG_Verified = .true.
      if (MaxDev > ThrDia*max(One,maxval(abs(DRef(:))))) then
        write(u6,*) 'Mk_IHI: sg_hdiag_fast differs from sg_hdiag_ref'
        call Abend()
      end if
      call mma_deallocate(DRef)
    end if
    call mma_deallocate(DAvg)
  end subroutine SGCI_Diag

  subroutine Verify_IHI()
    use rasscf_global, only: NSEL
    integer(kind=iwp), parameter :: MxRef = 500, MxPrt = 10
    integer(kind=iwp) :: i, iTyp, iCnf, k, nCnf, nC, iC, nSD, nSDErr, nSelL, nSelMin, nSelMean, nBoth, nPrt, iRk(3), iRkMin, &
                         iRkMax
    integer(kind=iwp), allocatable :: nCSFCnf(:), nOpCnf(:)
    real(kind=wp) :: Dev, MaxDev
    real(kind=wp), allocatable :: DSG(:), DRef(:), DEx(:), KeyL(:), KeyMin(:), KeyMean(:)
    logical(kind=iwp), allocatable :: SelL(:), SelMin(:), SelMean(:)

    ! exact diagonal in the SGUGA order
    call mma_allocate(DSG,nConf,Label='DSG')
    call SG_HDiag_G('FAST',iState,nConf,CI_St(iState)%STSym,TUVX,nTUVX,TU,nTU,DSG)
    if (nConf <= MxRef) then
      call mma_allocate(DRef,nConf,Label='DRef')
      call SG_HDiag_G('REF',iState,nConf,CI_St(iState)%STSym,TUVX,nTUVX,TU,nTU,DRef)
      MaxDev = Zero
      do i=1,nConf
        Dev = abs(DSG(i)-DRef(i))
        MaxDev = max(MaxDev,Dev)
        if (Dev > ThrDia*max(One,abs(DRef(i)))) then
          write(u6,*) 'Mk_IHI: sg_hdiag_fast differs from sg_hdiag_ref for CSF',i,DSG(i),DRef(i)
          call Abend()
        end if
      end do
      write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI: sg_hdiag_fast verified against sg_hdiag_ref, ',nConf,' CSFs, max. error ',MaxDev
      call mma_deallocate(DRef)
    end if

    ! to the Lucia CSF order (with HEXS/DEXS zero for the removed CSFs, as
    ! Lucia's)
    call HEXS_Proj(HEXS_Keep_SG,DSG)
    call mma_allocate(DEx,nConf,Label='DEx')
    call Diag_to_Lucia(DSG,DEx)
    call mma_deallocate(DSG)

    ! configuration layout (as Mk_Cnf_Layout) with the number of open shells
    ! (NDTFTP counts spin combinations, not determinants, for MS = 0)
    nCnf = sum(NCNFTP(1:NTYP,CI_St(iState)%STSym))
    call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
    call mma_allocate(nOpCnf,nCnf,Label='nOpCnf')
    k = 0
    do iTyp=1,NTYP
      do iCnf=1,NCNFTP(iTyp,CI_St(iState)%STSym)
        k = k+1
        nCSFCnf(k) = NCSFTP(iTyp)
        nOpCnf(k) = iTyp-1+MINOP
      end do
    end do

    ! keys per configuration: Lucia (the average, as PHPCSF uses it) and
    ! the exact lowest and mean CSF energies; exact check for single
    ! determinants: closed shells, or all open shells parallel (2S open
    ! shells, the high-spin determinant is the only one with MS = S)
    call mma_allocate(KeyL,nCnf,Label='KeyL')
    call mma_allocate(KeyMin,nCnf,Label='KeyMin')
    call mma_allocate(KeyMean,nCnf,Label='KeyMean')
    nSD = 0
    nSDErr = 0
    MaxDev = Zero
    iC = 0
    do k=1,nCnf
      nC = nCSFCnf(k)
      KeyL(k) = CSFDIA(iC+1)
      KeyMin(k) = minval(DEx(iC+1:iC+nC))
      KeyMean(k) = sum(DEx(iC+1:iC+nC))/real(nC,kind=wp)
      if (nOpCnf(k) == CI_St(iState)%iSpin-1) then
        nSD = nSD+1
        Dev = abs(CSFDIA(iC+1)-DEx(iC+1))
        MaxDev = max(MaxDev,Dev)
        if (Dev > ThrDia*max(One,abs(DEx(iC+1)))) then
          nSDErr = nSDErr+1
          if (nSDErr <= MxPrt) write(u6,*) 'Mk_IHI: single-determinant configuration',k,' CI solver:',CSFDIA(iC+1), &
                                           ' exact:',DEx(iC+1)
        end if
      end if
      iC = iC+nC
    end do
    if (nSDErr > 0) then
      write(u6,*) 'Mk_IHI: diagonal of single-determinant configurations differs from SGUGA'
      call Abend()
    end if
    write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI: diagonal verified against SGUGA for ',nSD,' single-determinant configurations, &
                                   &max. error ',MaxDev

    ! selection of PHPCSF with each key
    call mma_allocate(SelL,nCnf,Label='SelL')
    call mma_allocate(SelMin,nCnf,Label='SelMin')
    call mma_allocate(SelMean,nCnf,Label='SelMean')
    call Sel_PHP(KeyL,SelL)
    call Sel_PHP(KeyMin,SelMin)
    call Sel_PHP(KeyMean,SelMean)
    nSelL = count(SelL)
    nSelMin = count(SelMin)
    nSelMean = count(SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A)') 'Mk_IHI: explicit subspace (NSEL = ',NSEL,' of ',nConf,' CSFs, ',nCnf, &
                                           ' configurations) from the CI-solver diagonal: ',nSelL,' configurations'
    if (NSEL >= nConf) write(u6,'(1X,A)') '        (the whole space, the comparison is trivial)'
    nBoth = count(SelL .and. SelMin)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact lowest CSF energy per configuration: ',nSelMin,' configurations, ', &
                                      nBoth,' in common, ',count(SelL .neqv. SelMin),' differ'
    nBoth = count(SelL .and. SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact mean CSF energy per configuration:   ',nSelMean,' configurations, ', &
                                      nBoth,' in common, ',count(SelL .neqv. SelMean),' differ'
    ! The selection stops at the first configuration that does not fit into
    ! NSEL CSFs, so near the cutoff the sizes of the configurations decide;
    ! what matters is whether the differing configurations are ranked near
    ! the cutoff (rank = position in the order of increasing energy).
    if (any(SelL .neqv. SelMin) .or. any(SelL .neqv. SelMean)) then
      write(u6,'(1X,A)') '        differing configurations (selected by CI solver/exact lowest/exact mean, CSFs, '// &
                         'open shells, rank by CI-solver key, by exact lowest, by exact mean):'
      nPrt = 0
      iRkMin = huge(iRkMin)
      iRkMax = 0
      do k=1,nCnf
        if ((SelL(k) .eqv. SelMin(k)) .and. (SelL(k) .eqv. SelMean(k))) cycle
        iRk(1) = Cnf_Rank(KeyL,k)
        iRk(2) = Cnf_Rank(KeyMin,k)
        iRk(3) = Cnf_Rank(KeyMean,k)
        iRkMin = min(iRkMin,minval(iRk))
        iRkMax = max(iRkMax,maxval(iRk))
        nPrt = nPrt+1
        if (nPrt <= MxPrt) write(u6,'(9X,I8,3L2,2I6,3I8)') k,SelL(k),SelMin(k),SelMean(k),nCSFCnf(k),nOpCnf(k),iRk(:)
      end do
      write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        differing configurations are ranked ',iRkMin,' to ',iRkMax, &
                                        ' (',nSelL,' selected by the CI solver)'
    end if

    call mma_deallocate(SelMean)
    call mma_deallocate(SelMin)
    call mma_deallocate(SelL)
    call mma_deallocate(KeyMean)
    call mma_deallocate(KeyMin)
    call mma_deallocate(KeyL)
    call mma_deallocate(nOpCnf)
    call mma_deallocate(nCSFCnf)
    call mma_deallocate(DEx)

  end subroutine Verify_IHI

  ! A diagonal in the SGUGA CSF order to the Lucia order: a permutation,
  ! without the CSF phases (SG_ReOrd of the index vector: Idx(ICSF) = +-ISG).
  subroutine Diag_to_Lucia(DSG,DL)
    real(kind=wp), intent(in) :: DSG(nConf)
    real(kind=wp), intent(out) :: DL(nConf)
    integer(kind=iwp) :: j
    real(kind=wp), allocatable :: Idx(:), Map(:)

    call mma_allocate(Map,nConf,Label='Map')
    call mma_allocate(Idx,nConf,Label='Idx')
    do j=1,nConf
      Map(j) = real(j,kind=wp)
    end do
    call SG_ReOrd(iState,CI_St(iState)%STSym,1,nConf,Map,Idx)
    do j=1,nConf
      DL(j) = DSG(nint(abs(Idx(j))))
    end do
    call mma_deallocate(Idx)
    call mma_deallocate(Map)

  end subroutine Diag_to_Lucia

  ! Exact diagonal <I|H|I> in the Lucia CSF order (sg_hdiag_fast).
  subroutine Exact_Diag(DL)
    real(kind=wp), intent(out) :: DL(nConf)
    real(kind=wp), allocatable :: DSG(:)

    call mma_allocate(DSG,nConf,Label='DSG')
    call SG_HDiag_G('FAST',iState,nConf,CI_St(iState)%STSym,TUVX,nTUVX,TU,nTU,DSG)
    call HEXS_Proj(HEXS_Keep_SG,DSG)
    call Diag_to_Lucia(DSG,DL)
    call mma_deallocate(DSG)

  end subroutine Exact_Diag

  ! Position of Key(k) in the order of increasing Key (ties by index).
  function Cnf_Rank(Key,k) result(iRank)
    real(kind=wp), intent(in) :: Key(:)
    integer(kind=iwp), intent(in) :: k
    integer(kind=iwp) :: iRank, l

    iRank = 1
    do l=1,size(Key)
      if ((Key(l) < Key(k)) .or. ((Key(l) == Key(k)) .and. (l < k))) iRank = iRank+1
    end do

  end function Cnf_Rank

  ! Configurations selected by PHPCSF for the explicit subspace with the
  ! configuration energies Key: the lowest first (ties within Acc in the
  ! order of the configurations), until the next one does not fit into NSEL
  ! CSFs.
  subroutine Sel_PHP(Key,Sel)
    use rasscf_global, only: NSEL
    real(kind=wp), intent(in) :: Key(:)
    logical(kind=iwp), intent(out) :: Sel(size(Key))
    integer(kind=iwp) :: k, kMin, nP, nPCnf, nCMin
    real(kind=wp) :: XMax, XMin
    real(kind=wp), allocatable :: Wrk(:)
    real(kind=wp), parameter :: Acc = 1.0e-13_wp
    integer(kind=iwp), allocatable :: nCSFCnf(:)

    call Mk_Cnf_Layout_Lucia(CI_St(iState)%STSym,nCSFCnf)
    call mma_allocate(Wrk,size(Key),Label='Wrk')
    Wrk(:) = Key(:)
    XMax = maxval(abs(Key))
    Sel(:) = .false.
    nP = 0
    nPCnf = 0
    do
      XMin = XMax+One
      kMin = 0
      nCMin = 0
      do k=1,size(Key)
        if (Wrk(k)+Acc < XMin) then
          XMin = Wrk(k)
          kMin = k
          nCMin = nCSFCnf(k)
        end if
      end do
      if ((kMin == 0) .or. (nP+nCMin > NSEL)) exit
      Sel(kMin) = .true.
      nP = nP+nCMin
      nPCnf = nPCnf+1
      Wrk(kMin) = XMax+One
      if (nPCnf >= size(Key)) exit
    end do
    call mma_deallocate(Wrk)
    call mma_deallocate(nCSFCnf)

  end subroutine Sel_PHP

End subroutine Mk_IHI

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Layout of the CSF basis of symmetry IREFSM by configuration:
!> nCSFCnf(k) is the number of CSFs of the k-th configuration, the
!> configurations being in the order of the CSF basis (the CSFs of one
!> configuration are consecutive). Allocated here, to be freed by the caller
!> with mma_deallocate. (Lucia back end: configuration types from lucia_data.)
!> Layout of the CSF basis of symmetry IREFSM by configuration: nCSFCnf(k)
!> CSFs in configuration k, and optionally their indices in the CSF basis,
!> iCSFCnf(iOff+1:iOff+nCSFCnf(k)) (configuration by configuration). With
!> SGCI the SGUGA configurations and split-graph CSF indices (not
!> contiguous), otherwise those of the CI back end (contiguous).
subroutine Mk_Cnf_Layout(IREFSM,nCSFCnf,iCSFCnf)
integer(kind=iwp), intent(in) :: IREFSM
integer(kind=iwp), allocatable, intent(out) :: nCSFCnf(:)
integer(kind=iwp), allocatable, intent(out), optional :: iCSFCnf(:)
integer(kind=iwp) :: k, I

if (SG_CI()) then
  call SG_Cnf_Kept(IREFSM)
  call mma_allocate(nCSFCnf,Cnf_SG%nCnf,Label='nCSFCnf')
  do k=1,Cnf_SG%nCnf
    nCSFCnf(k) = Cnf_SG%iOff(k+1)-Cnf_SG%iOff(k)
  end do
  if (present(iCSFCnf)) then
    call mma_allocate(iCSFCnf,Cnf_SG%iOff(Cnf_SG%nCnf+1),Label='iCSFCnf')
    iCSFCnf(:) = Cnf_SG%iCSF(1:Cnf_SG%iOff(Cnf_SG%nCnf+1))
  end if
else if (DoFaro) then
  ! Faroald (citrans): the CSFs of a configuration are contiguous
  call Mk_Cnf_Layout_Faroald(nCSFCnf)
  if (present(iCSFCnf)) then
    call mma_allocate(iCSFCnf,sum(nCSFCnf),Label='iCSFCnf')
    do I=1,size(iCSFCnf)
      iCSFCnf(I) = I
    end do
  end if
else
  call Mk_Cnf_Layout_Lucia(IREFSM,nCSFCnf)
  if (present(iCSFCnf)) then
    call mma_allocate(iCSFCnf,sum(nCSFCnf),Label='iCSFCnf')
    do I=1,size(iCSFCnf)
      iCSFCnf(I) = I
    end do
  end if
end if

end subroutine Mk_Cnf_Layout

!> The SGUGA configuration list of symmetry LSM, built once and kept
!> (Cnf_SG) for SGCI.
subroutine SG_Cnf_Kept(LSM)
use sguga, only: SG_Mk_Cnf_List
integer(kind=iwp), intent(in) :: LSM
! the CI state of RASSCF (SGS, CIS, EXS)
integer(kind=iwp), parameter :: iState = 1
if ((Cnf_SG%nCnf == 0) .or. (Cnf_SG%LSM /= LSM)) call SG_Mk_Cnf_List(iState,LSM,Cnf_SG)
end subroutine SG_Cnf_Kept

!> The configurations of Faroald (citrans order) and their numbers of CSFs.
subroutine Mk_Cnf_Layout_Faroald(nCSFCnf)
use citrans, only: citrans_conf_dets
integer(kind=iwp), allocatable, intent(out) :: nCSFCnf(:)
integer(kind=iwp) :: iPass, k, k1, nc, nd, nCnf
integer(kind=iwp), allocatable :: da(:), db(:)
real(kind=wp), allocatable :: cf(:,:)
do iPass=1,2
  nCnf = 0
  k = 1
  do while (k <= nConf)
    call citrans_conf_dets(k,k1,nc,nd,da,db,cf)
    nCnf = nCnf+1
    if (iPass == 2) nCSFCnf(nCnf) = nc
    k = k1+nc
    call mma_deallocate(da)
    call mma_deallocate(db)
    call mma_deallocate(cf)
  end do
  if (iPass == 1) call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
end do
end subroutine Mk_Cnf_Layout_Faroald

subroutine Mk_Cnf_Layout_Lucia(IREFSM,nCSFCnf)
use lucia_data, only: NCNFTP, NCSFTP, NTYP

integer(kind=iwp), intent(in) :: IREFSM
integer(kind=iwp), allocatable, intent(out) :: nCSFCnf(:)
integer(kind=iwp) :: iTyp, iCnf, k

call mma_allocate(nCSFCnf,sum(NCNFTP(1:NTYP,IREFSM)),Label='nCSFCnf')
k = 0
do iTyp=1,NTYP
  do iCnf=1,NCNFTP(iTyp,IREFSM)
    k = k+1
    nCSFCnf(k) = NCSFTP(iTyp)
  end do
end do

end subroutine Mk_Cnf_Layout_Lucia

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Explicit CI Hamiltonian <I|H|J> over the CSFs of the configurations
!> IPCNF(1:NPCNF) of symmetry IREFSM (configuration indices as in
!> Mk_Cnf_Layout), packed as a lower triangle: the CSFs are ordered
!> configuration by configuration in the order given, the CSFs of each
!> configuration consecutively. TU and TUVX are the one- and two-electron
!> integrals in triangular packing (as for Mk_IHI); ECORE is added on the
!> diagonal. The Hamiltonian is defined by TU and TUVX alone, as for the
!> sigma vector (Mk_H_Psi): no scaling of the exchange (ExFac), which would
!> only be consistent if built into the integrals by the caller.
!>
!> Back ends:
!>  - Lucia: CNHCN for each pair of configurations (determinant block,
!>    transformed with DTOC).
!>  - Faroald (DoFaro): for each pair of configurations the determinant
!>    block (Slater-Condon rules, Faroald_HDet) transformed with the spin
!>    tables of citrans, C_L^T H C_R.
!>  - SGUGA (RASSCF keyword SGEX; no GAS, single RAS partition): <I|H|J>
!>    for each pair of CSFs from the occupation difference (SG_MatEl_H),
!>    pairs of configurations more than a double excitation apart skipped.
!> With MOLCAS_CI_VERIFY=SIGMA, once per run, the block in use is compared
!> with the other back ends (the Faroald block also with columns H e_J of
!> the Faroald sigma vector) and checked against Mk_H_Psi; the CPU time of
!> each back end is reported.
subroutine Mk_IHJ(PHP,nPHP,IPCNF,NPCNF,IREFSM,ECORE,TU,nTU,TUVX,nTUVX,NTEST)
use Index_Functions, only: nTri_Elem
use Lucia_data, only: Conf_Occ, DFTP, DTOC, IREOTS, Sigma_on_disk
use csfbas, only: NAEL, NBEL
use lucia_data, only: NCSFTP, NDTFTP, NTYP
use Constants, only: One

integer(kind=iwp), intent(in) :: nPHP, NPCNF, IPCNF(NPCNF), IREFSM, nTU, nTUVX
real(kind=wp), intent(out) :: PHP(nPHP)
real(kind=wp), intent(in) :: ECORE, TU(nTU), TUVX(nTUVX)
integer(kind=iwp), intent(inout) :: NTEST
integer(kind=iwp) :: NPCSF
integer(kind=iwp), allocatable :: IPCSF(:), nSelCnf(:)
real(kind=wp) :: T0, T1, W
real(kind=wp), allocatable :: PHP_Ref(:)
real(kind=wp), parameter :: Thr = 1.0e-10_wp
! the CI state of RASSCF (SGS, CIS, EXS)
integer(kind=iwp), parameter :: iState = 1
logical(kind=iwp) :: SG_OK, UseSG, Verify
character(len=8) :: Lab

! CSF indices (in the CSF basis of symmetry IREFSM) of the selection.
call Sel_CSFs()
NPCSF = size(IPCSF)
if (nTri_Elem(NPCSF) > nPHP) then
  write(u6,*) 'Mk_IHJ: PHP too small for the selected CSFs:',nTri_Elem(NPCSF),nPHP
  call Abend()
end if

! SGCI: SGUGA configurations and split-graph CSFs; under SIGMA the block is
! checked once against exact columns H e_J of Mk_H_Psi (itself compared
! with Lucia there).
if (SG_CI()) then
  call SG_Cnf_Kept(IREFSM)
  call SG_Block_G(Cnf_SG,NPCNF,IPCNF,PHP)
  if (CI_Verify('SIGMA') .and. (.not. IHJ_Verified)) then
    IHJ_Verified = .true.
    call Verify_IHJ()
  end if
  call mma_deallocate(IPCSF)
  call mma_deallocate(nSelCnf)
  return
end if

! SGUGA is available without GAS and for a single RAS partition (SGS is not
! set up with GAS, so test iDoGAS first).
SG_OK = SG_Avail()
if (DoSGEx .and. (.not. SG_OK) .and. (.not. SGEx_Warned)) then
  write(u6,'(1X,A)') 'Mk_IHJ: SGEX ignored, SGUGA is not available with GAS; using the CI back end'
  SGEx_Warned = .true.
end if
! FARO: the explicit block is always that of Faroald, in its own basis
! (SGEX belongs to the Lucia path)
if (DoSGEx .and. DoFaro .and. (.not. SGEx_Warned)) then
  write(u6,'(1X,A)') 'Mk_IHJ: SGEX ignored with FARO; the explicit block is computed by Faroald'
  SGEx_Warned = .true.
end if
UseSG = DoSGEx .and. SG_OK .and. (.not. DoFaro)
Verify = CI_Verify('SIGMA') .and. (.not. IHJ_Verified)
! the walks of the Lucia CSFs, once per run
if (Verify .and. SG_OK .and. (.not. DoFaro)) call SG_Verify_Lucia_Walk(iState,IREFSM,nConf)

call CWTime(T0,W)
if (UseSG) then
  Lab = 'SGUGA'
  call IHJ_SGUGA(PHP)
else if (DoFaro) then
  Lab = 'Faroald'
  call IHJ_Faroald(PHP)
else
  Lab = 'Lucia'
  call IHJ_Lucia(PHP)
end if
call CWTime(T1,W)

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare the block, once per run,
! with the other back ends and with exact columns H e_J from Mk_H_Psi.
if (Verify) then
  write(u6,'(1X,3A,I0,A,F10.3,A)') 'Mk_IHJ: ',trim(Lab),' block, ',NPCSF,' CSFs, CPU ',T1-T0,' s'
  call mma_allocate(PHP_Ref,nTri_Elem(NPCSF),Label='PHP_Ref')
  ! (with FARO the block is in the basis of Faroald: the references of
  ! Lucia and SGUGA are their columns H e_J in that basis, Ref_Sigma)
  if (UseSG .and. (.not. DoFaro)) call Ref_Block('Lucia')
  if (DoFaro) then
    call Ref_Block('Faroald H e_J')
    call Ref_Block('Lucia H e_J')
    if (SG_OK) call Ref_Block('SGUGA H e_J')
  end if
  if (SG_OK) then
    if ((.not. UseSG) .and. (.not. DoFaro)) call Ref_Block('SGUGA')
    call Verify_SG_Explicit()
  else
    write(u6,'(1X,A)') 'Mk_IHJ: SGUGA block skipped (not available with GAS)'
  end if
  call mma_deallocate(PHP_Ref)
  IHJ_Verified = .true.
  call Verify_IHJ()
end if

call mma_deallocate(IPCSF)
call mma_deallocate(nSelCnf)

contains

  ! SGUGA explicit block over the configurations IPC(1:NPC) of Cnf; for GAS
  ! the matrix elements in the relaxed space
  subroutine SG_Block_G(Cnf,NPC,IPC,HB)
    type(SGCnfList), intent(in) :: Cnf
    integer(kind=iwp), intent(in) :: NPC, IPC(NPC)
    real(kind=wp), intent(out) :: HB(*)
    if (CI_St(iState)%iRel > 0) then
      call SG_Rel_Map(iState,Cnf%LSM)
      call SG_Explicit_Block(iState,Cnf,NPC,IPC,ECORE,TU,nTU,TUVX,nTUVX,HB,iRel=CI_St(iState)%iRel,RelMap=Rel_Map(Cnf%LSM)%M)
    else
      call SG_Explicit_Block(iState,Cnf,NPC,IPC,ECORE,TU,nTU,TUVX,nTUVX,HB)
    end if
  end subroutine SG_Block_G

  ! Reference block of back end LabRef in PHP_Ref, compared with PHP.
  subroutine Ref_Block(LabRef)
    character(len=*), intent(in) :: LabRef
    real(kind=wp) :: T2, T3, W2

    call CWTime(T2,W2)
    select case (LabRef)
      case ('Lucia')
        call IHJ_Lucia(PHP_Ref)
      case ('Faroald')
        call IHJ_Faroald(PHP_Ref)
      case ('Faroald H e_J')
        call IHJ_Faroald_Col(PHP_Ref)
      case ('Lucia H e_J')
        call IHJ_Ref_Col('Lucia',PHP_Ref)
      case ('SGUGA H e_J')
        call IHJ_Ref_Col('SGUGA',PHP_Ref)
      case ('SGUGA')
        call IHJ_SGUGA(PHP_Ref)
    end select
    call CWTime(T3,W2)
    write(u6,'(1X,3A,F10.3,A)') 'Mk_IHJ: ',LabRef,' block, CPU ',T3-T2,' s'
    call Compare_IHJ(trim(Lab),LabRef,PHP_Ref)

  end subroutine Ref_Block

  ! IPCSF: the CSFs of the selected configurations, in the order of PHP;
  ! nSelCnf: the number of CSFs of each selected configuration.
  subroutine Sel_CSFs()
    integer(kind=iwp) :: iCnf, iOff, kCnf, I, ICSF
    integer(kind=iwp), allocatable :: nCSFCnf(:), iCSF1(:), iCSFCnf(:)

    call Mk_Cnf_Layout(IREFSM,nCSFCnf,iCSFCnf)
    ! iCSF1(k): position of the first CSF of configuration k in iCSFCnf
    call mma_allocate(iCSF1,size(nCSFCnf),Label='iCSF1')
    iOff = 1
    do iCnf=1,size(nCSFCnf)
      iCSF1(iCnf) = iOff
      iOff = iOff+nCSFCnf(iCnf)
    end do
    call mma_allocate(IPCSF,sum(nCSFCnf(IPCNF(1:NPCNF))),Label='IPCSF')
    call mma_allocate(nSelCnf,NPCNF,Label='nSelCnf')
    ICSF = 0
    do kCnf=1,NPCNF
      nSelCnf(kCnf) = nCSFCnf(IPCNF(kCnf))
      do I=1,nCSFCnf(IPCNF(kCnf))
        ICSF = ICSF+1
        IPCSF(ICSF) = iCSFCnf(iCSF1(IPCNF(kCnf))+I-1)
      end do
    end do
    call mma_deallocate(iCSFCnf)
    call mma_deallocate(iCSF1)
    call mma_deallocate(nCSFCnf)

  end subroutine Sel_CSFs

  ! Lucia back end: CNHCN for each pair of selected configurations.
  subroutine IHJ_Lucia(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: NEL, MXDT, MXCS, nScr, ICNL, ICNR, ILTYP, IRTYP, NCSFL, NCSFR, IILB, IIRB, IIL, IIR, IIRMAX, DL, SL, DR, &
                         SR
    integer(kind=iwp), allocatable :: ICONFL(:), ICONFR(:)
    real(kind=wp), allocatable :: HONE(:,:), CNHCNM(:), Scr(:)

    NEL = NAEL+NBEL
    MXDT = maxval(NDTFTP(1:NTYP))
    MXCS = maxval(NCSFTP(1:NTYP))

    ! One-body matrix in square storage, as CNHCN wants it.
    call mma_allocate(HONE,NAC,NAC,Label='HONE')
    call SQUARE(TU,HONE,NAC,1,NAC)

    ! Scratch for CNHCN: determinant strings and signs of both configurations,
    ! CNFSTR work space, the determinant block, DIHDJ_MOLCAS work space and the
    ! half-transformed block.
    nScr = 2*MXDT*(NEL+1)+NEL*(MXDT+1)+MXDT*MXDT+4*NAC+MXDT+MXCS*MXDT+64
    call mma_allocate(Scr,nScr,Label='IHJ_Scr')
    call mma_allocate(CNHCNM,MXCS*MXCS,Label='CNHCNM')
    call mma_allocate(ICONFL,NEL,Label='ICONFL')
    call mma_allocate(ICONFR,NEL,Label='ICONFR')

    IILB = 1
    do ICNL=1,NPCNF
      call GETCNF(ICONFL,ILTYP,IPCNF(ICNL),Conf_Occ(IREFSM)%A,IREFSM,NEL)
      NCSFL = NCSFTP(ILTYP)
      call Lucia_Occ(ICONFL,ILTYP,DL,SL)
      IIRB = 1
      do ICNR=1,ICNL
        call GETCNF(ICONFR,IRTYP,IPCNF(ICNR),Conf_Occ(IREFSM)%A,IREFSM,NEL)
        NCSFR = NCSFTP(IRTYP)
        call Lucia_Occ(ICONFR,IRTYP,DR,SR)
        if (Cnf_Exc(DL,SL,DR,SR) > 2) then
          CNHCNM(1:NCSFL*NCSFR) = Zero
        else
          call CNHCN(ICONFL,ILTYP,ICONFR,IRTYP,CNHCNM,Scr,NAEL,NBEL,ECORE,HONE,DFTP,DTOC,NAC,TUVX,NTEST,One,IREOTS)
        end if
        do IIL=1,NCSFL
          IIRMAX = merge(IIL,NCSFR,IILB == IIRB)
          do IIR=1,IIRMAX
            H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = CNHCNM((IIR-1)*NCSFL+IIL)
          end do
        end do
        IIRB = IIRB+NCSFR
      end do
      IILB = IILB+NCSFL
    end do

    call mma_deallocate(ICONFR)
    call mma_deallocate(ICONFL)
    call mma_deallocate(CNHCNM)
    call mma_deallocate(Scr)
    call mma_deallocate(HONE)

  end subroutine IHJ_Lucia

  ! Doubly (D) and singly (S) occupied orbitals, as bit strings (orbital t in
  ! bit t-1), of the Lucia configuration ICONF of type ITYP (GETCNF form:
  ! the closed orbitals first, then the open ones).
  subroutine Lucia_Occ(ICONF,ITYP,D,S)
    use lucia_data, only: MINOP
    integer(kind=iwp), intent(in) :: ICONF(*), ITYP
    integer(kind=iwp), intent(out) :: D, S
    integer(kind=iwp) :: IOP, ICL, I

    IOP = ITYP-1+MINOP
    ICL = (NAEL+NBEL-IOP)/2
    D = 0
    S = 0
    do I=1,ICL
      D = ibset(D,ICONF(I)-1)
    end do
    do I=ICL+1,ICL+IOP
      S = ibset(S,ICONF(I)-1)
    end do

  end subroutine Lucia_Occ

  ! Faroald back end: for each pair of selected configurations L, R the
  ! determinant block <D_L|H|D_R> (Faroald_HDet) is transformed to the CSFs
  ! with the spin tables of citrans, C_L^T H C_R; pairs of configurations
  ! more than a double excitation apart are skipped (Cnf_Exc). The CI
  ! vectors (and the selected CSFs) are in the configuration order of
  ! citrans (FARO).
  subroutine IHJ_Faroald(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, kCnf, ICNL, ICNR, IILB, IIRB, IIL, IIR, IIRMAX, k1, nc, nd, nDTot, nCTot, nL, nR, ndL, ndR
    integer(kind=iwp), allocatable :: kSel(:), iSel(:), iDOff(:), iCOff(:), nDCnf(:), DA(:), DB(:), da1(:), db1(:), &
                                      DOcc(:), SOcc(:)
    real(kind=wp), allocatable :: C(:), cf(:,:), HDet(:,:), T(:,:), Blk(:,:)

    call Faroald_Ints(TU,TUVX)

    ! selected CSF J <-> citrans CSF kSel(J), with phase iSel(J)
    call mma_allocate(kSel,NPCSF,Label='kSel')
    call mma_allocate(iSel,NPCSF,Label='iSel')
    ! (the selected CSFs are citrans CSFs: FARO works in that basis)
    kSel(:) = IPCSF(:)
    iSel(:) = 1

    ! Determinants and CSF coefficients of the selected configurations.
    call mma_allocate(nDCnf,NPCNF,Label='nDCnf')
    call mma_allocate(iDOff,NPCNF,Label='iDOff')
    call mma_allocate(iCOff,NPCNF,Label='iCOff')
    call mma_allocate(DOcc,NPCNF,Label='DOcc')
    call mma_allocate(SOcc,NPCNF,Label='SOcc')
    nDTot = 0
    nCTot = 0
    IILB = 0
    do kCnf=1,NPCNF
      call citrans_conf_dets(kSel(IILB+1),k1,nc,nd,da1,db1,cf)
      if (nc /= nSelCnf(kCnf)) then
        write(u6,*) 'Mk_IHJ: different number of CSFs in Lucia and citrans:',nSelCnf(kCnf),nc
        call Abend()
      end if
      nDCnf(kCnf) = nd
      ! doubly and singly occupied orbitals of the configuration
      DOcc(kCnf) = iand(da1(1),db1(1))
      SOcc(kCnf) = ieor(da1(1),db1(1))
      iDOff(kCnf) = nDTot
      iCOff(kCnf) = nCTot
      nDTot = nDTot+nd
      nCTot = nCTot+nd*nc
      IILB = IILB+nc
      call mma_deallocate(da1)
      call mma_deallocate(db1)
      call mma_deallocate(cf)
    end do
    call mma_allocate(DA,nDTot,Label='DA')
    call mma_allocate(DB,nDTot,Label='DB')
    call mma_allocate(C,nCTot,Label='C')
    IILB = 0
    do kCnf=1,NPCNF
      call citrans_conf_dets(kSel(IILB+1),k1,nc,nd,da1,db1,cf)
      DA(iDOff(kCnf)+1:iDOff(kCnf)+nd) = da1(:)
      DB(iDOff(kCnf)+1:iDOff(kCnf)+nd) = db1(:)
      do I=1,nc
        J = kSel(IILB+I)-k1+1
        if ((J < 1) .or. (J > nc)) then
          write(u6,*) 'Mk_IHJ: CSF of the selection not in its configuration:',IPCSF(IILB+I)
          call Abend()
        end if
        C(iCOff(kCnf)+(I-1)*nd+1:iCOff(kCnf)+I*nd) = iSel(IILB+I)*cf(:,J)
      end do
      IILB = IILB+nc
      call mma_deallocate(da1)
      call mma_deallocate(db1)
      call mma_deallocate(cf)
    end do
    call mma_deallocate(iSel)
    call mma_deallocate(kSel)

    ! Blocks C_L^T <D_L|H|D_R> C_R.
    nd = maxval(nDCnf)
    nc = maxval(nSelCnf)
    call mma_allocate(HDet,nd,nd,Label='HDet')
    call mma_allocate(T,nd,nc,Label='T')
    call mma_allocate(Blk,nc,nc,Label='Blk')
    IILB = 1
    do ICNL=1,NPCNF
      nL = nSelCnf(ICNL)
      ndL = nDCnf(ICNL)
      IIRB = 1
      do ICNR=1,ICNL
        nR = nSelCnf(ICNR)
        ndR = nDCnf(ICNR)
        if (Cnf_Exc(DOcc(ICNL),SOcc(ICNL),DOcc(ICNR),SOcc(ICNR)) > 2) then
          Blk(1:nL,1:nR) = Zero
        else
          do J=1,ndR
            do I=1,ndL
              HDet(I,J) = Faroald_HDet(DA(iDOff(ICNL)+I),DB(iDOff(ICNL)+I),DA(iDOff(ICNR)+J),DB(iDOff(ICNR)+J))
            end do
          end do
          call DGEMM_('N','N',ndL,nR,ndR,One,HDet,nd,C(iCOff(ICNR)+1:iCOff(ICNR)+ndR*nR),ndR,Zero,T,nd)
          call DGEMM_('T','N',nL,nR,ndL,One,C(iCOff(ICNL)+1:iCOff(ICNL)+ndL*nL),ndL,T,nd,Zero,Blk,nc)
        end if
        do IIL=1,nL
          IIRMAX = merge(IIL,nR,IILB == IIRB)
          do IIR=1,IIRMAX
            H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = Blk(IIL,IIR)
          end do
        end do
        IIRB = IIRB+nR
      end do
      IILB = IILB+nL
    end do
    do I=1,NPCSF
      H(nTri_Elem(I)) = H(nTri_Elem(I))+ECORE
    end do

    call mma_deallocate(Blk)
    call mma_deallocate(T)
    call mma_deallocate(HDet)
    call mma_deallocate(C)
    call mma_deallocate(DB)
    call mma_deallocate(DA)
    call mma_deallocate(SOcc)
    call mma_deallocate(DOcc)
    call mma_deallocate(iCOff)
    call mma_deallocate(iDOff)
    call mma_deallocate(nDCnf)

  end subroutine IHJ_Faroald

  ! SGUGA back end: <I|H|J> = P_I P_J <I_G|H|J_G> + delta_IJ ECORE, with
  ! I_G the GUGA CSF of the Lucia CSF I and P_I the phase of SG_ReOrd
  ! (mode 0: c_G(I_G) = P_I c(I)); the block <I_G|H|J_G> from SG_Block_H
  ! (pairs of configurations more than a double excitation apart skipped).
  subroutine IHJ_SGUGA(H)
    use Molcas, only: MxAct
    use sguga, only: SG_Block_H, SG_MatEl_Init, SG_NUM
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, iS
    integer(kind=iwp) :: IWALK(MxAct)
    integer(kind=iwp), allocatable :: iG(:), iPG(:)

    ! GAS: in the relaxed space (the CSFs by their walks)
    iS = iState
    if (CI_St(iState)%iRel > 0) iS = CI_St(iState)%iRel
    call SG_MatEl_Init(iS,IREFSM)
    call mma_allocate(iG,NPCSF,Label='iG')
    call mma_allocate(iPG,NPCSF,Label='iPG')
    do J=1,NPCSF
      call SG_Lucia_Walk(iState,IREFSM,IPCSF(J),IWALK,iPG(J))
      iG(J) = SG_NUM(SGS(iS),EXS(iS),IWALK)
    end do
    call SG_Block_H(IREFSM,NPCNF,nSelCnf,NPCSF,iG,TU,nTU,TUVX,nTUVX,H)
    do I=1,NPCSF
      do J=1,I
        H(nTri_Elem(I-1)+J) = iPG(I)*iPG(J)*H(nTri_Elem(I-1)+J)
      end do
      H(nTri_Elem(I)) = H(nTri_Elem(I))+ECORE
    end do
    call mma_deallocate(iPG)
    call mma_deallocate(iG)

  end subroutine IHJ_SGUGA

  ! CI development (SIGMA): SGUGA's own configuration list and explicit
  ! subspace (SG_Mk_Cnf_List, SG_Select_Explicit, SG_Explicit_Block) against
  ! the CI back end: (1) each configuration of the back end is exactly one
  ! SGUGA configuration (same CSFs), (2) the two selections are compared
  ! (the keys differ: the configuration average of the back end, the mean
  ! exact CSF diagonal of SGUGA), (3) the SGUGA block over the configurations
  ! selected by the back end has the eigenvalues of the block in use.
  subroutine Verify_SG_Explicit()
    use Constants, only: Two
    use sguga, only: SGCnfList, SG_Mk_Cnf_List, SG_Free_Cnf_List, SG_Select_Explicit, SG_Cnf_Avg_Diag
    type(SGCnfList) :: Cnf
    integer(kind=iwp) :: I, k, kG, iOffL, MxDim, nGSel, nGCSF, nBoth, nLSel
    integer(kind=iwp), allocatable :: nCSFL(:), LtoG(:), CnfOfG(:), LtoGCnf(:), GSel(:), GCSF(:), GtoL(:)
    logical(kind=iwp), allocatable :: SelG(:), SelL(:)
    real(kind=wp) :: Dev
    real(kind=wp), allocatable :: Idx(:), Map(:), DSG(:), HG(:), EG(:), EL(:)

    call SG_Mk_Cnf_List(iState,IREFSM,Cnf)
    ! back-end CSF -> split-graph CSF, and the SGUGA configuration of each
    call mma_allocate(Idx,nConf,Label='Idx')
    call mma_allocate(Map,nConf,Label='Map')
    do I=1,nConf
      Idx(I) = real(I,kind=wp)
    end do
    call CI_Basis_To_File(iState,IREFSM,nConf,Idx,Map)
    call mma_allocate(LtoG,nConf,Label='LtoG')
    do I=1,nConf
      LtoG(nint(abs(Map(I)))) = I
    end do
    call mma_deallocate(Map)
    call mma_deallocate(Idx)
    call mma_allocate(CnfOfG,nConf,Label='CnfOfG')
    do kG=1,Cnf%nCnf
      CnfOfG(Cnf%iCSF(Cnf%iOff(kG)+1:Cnf%iOff(kG+1))) = kG
    end do
    ! (1) configurations of the back end -> SGUGA configurations, a bijection
    call Mk_Cnf_Layout(IREFSM,nCSFL)
    if (size(nCSFL) /= Cnf%nCnf) then
      write(u6,*) 'Mk_IHJ: number of configurations, CI back end and SGUGA:',size(nCSFL),Cnf%nCnf
      call Abend()
    end if
    call mma_allocate(LtoGCnf,size(nCSFL),Label='LtoGCnf')
    call mma_allocate(SelG,Cnf%nCnf,Label='SelG')
    SelG(:) = .false.
    iOffL = 0
    do k=1,size(nCSFL)
      kG = CnfOfG(LtoG(iOffL+1))
      if (any(CnfOfG(LtoG(iOffL+1:iOffL+nCSFL(k))) /= kG) .or. (nCSFL(k) /= Cnf%iOff(kG+1)-Cnf%iOff(kG)) .or. SelG(kG)) then
        write(u6,*) 'Mk_IHJ: configuration',k,' of the CI back end is not one SGUGA configuration'
        call Abend()
      end if
      SelG(kG) = .true.
      LtoGCnf(k) = kG
      iOffL = iOffL+nCSFL(k)
    end do
    write(u6,'(1X,A,I0,A)') 'Mk_IHJ: SGUGA configuration list verified, ',Cnf%nCnf,' configurations as in the CI back end'
    ! (2) the SGUGA selection with its configuration-average diagonal (the
    ! key of SGCI and of Lucia: differences come from ties only)
    MxDim = nint((sqrt(real(8*nPHP+1,kind=wp))-One)/Two)
    call mma_allocate(DSG,nConf,Label='DSG')
    call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,DSG)
    call HEXS_Proj(HEXS_Keep_SG,DSG)
    call mma_allocate(GSel,Cnf%nCnf,Label='GSel')
    call mma_allocate(GCSF,max(1,MxDim),Label='GCSF')
    call SG_Select_Explicit(Cnf,DSG,MxDim,nGSel,GSel,nGCSF,GCSF)
    call mma_deallocate(DSG)
    call mma_allocate(SelL,Cnf%nCnf,Label='SelL')
    SelG(:) = .false.
    SelL(:) = .false.
    SelG(GSel(1:nGSel)) = .true.
    SelL(LtoGCnf(IPCNF(1:NPCNF))) = .true.
    nBoth = count(SelG .and. SelL)
    nLSel = count(SelL)
    write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'Mk_IHJ: SGUGA selection (configuration average): ',nGSel, &
          ' configurations, ',nGCSF,' CSFs; CI back end: ',nLSel,' configurations; ',nBoth,' in common, ', &
          count(SelG .neqv. SelL),' differ'
    ! (3) SGUGA block over the configurations selected by the back end
    call mma_allocate(GtoL,max(1,NPCNF),Label='GtoL')
    GtoL(1:NPCNF) = LtoGCnf(IPCNF(1:NPCNF))
    call mma_allocate(HG,nTri_Elem(NPCSF),Label='HG')
    call SG_Block_G(Cnf,NPCNF,GtoL,HG)
    call mma_allocate(EG,NPCSF,Label='EG')
    call mma_allocate(EL,NPCSF,Label='EL')
    call Block_Eigenvalues(HG,NPCSF,EG)
    HG(:) = PHP(1:nTri_Elem(NPCSF))
    call Block_Eigenvalues(HG,NPCSF,EL)
    Dev = maxval(abs(EG(:)-EL(:)))
    write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHJ: SG_Explicit_Block, eigenvalues verified against the block in use, ',NPCSF, &
                                   ' CSFs, max. error ',Dev
    if (Dev > Thr*max(One,maxval(abs(EL(:))))) then
      write(u6,*) 'Mk_IHJ: eigenvalues of the SGUGA explicit block differ'
      call Abend()
    end if
    call mma_deallocate(EL)
    call mma_deallocate(EG)
    call mma_deallocate(HG)
    call mma_deallocate(GtoL)
    call mma_deallocate(SelL)
    call mma_deallocate(GCSF)
    call mma_deallocate(GSel)
    call mma_deallocate(SelG)
    call mma_deallocate(LtoGCnf)
    call mma_deallocate(nCSFL)
    call mma_deallocate(CnfOfG)
    call mma_deallocate(LtoG)
    call SG_Free_Cnf_List(Cnf)

  end subroutine Verify_SG_Explicit

  ! Eigenvalues (increasing) of a lower-triangle packed symmetric matrix;
  ! the matrix is overwritten.
  subroutine Block_Eigenvalues(A,n,E)
    integer(kind=iwp), intent(in) :: n
    real(kind=wp), intent(inout) :: A(nTri_Elem(n))
    real(kind=wp), intent(out) :: E(n)
    integer(kind=iwp) :: I
    real(kind=wp), allocatable :: V(:,:)
    call mma_allocate(V,n,n,Label='V')
    V(:,:) = Zero
    do I=1,n
      V(I,I) = One
    end do
    call NIdiag_new(A,V,n,n)
    call JACORD(A,V,n,n)
    do I=1,n
      E(I) = A(nTri_Elem(I))
    end do
    call mma_deallocate(V)
  end subroutine Block_Eigenvalues

  ! Faroald reference (CI development): H(I,J) = (H e_J)(I) from the Faroald
  ! sigma vector of the unit vector e_J, J over the selected CSFs.
  subroutine IHJ_Faroald_Col(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, nTmp
    real(kind=wp), allocatable :: EJ(:), HEJ(:), cTmp(:), sTmp(:)

    nTmp = Faro_nTmp()
    call mma_allocate(cTmp,nTmp,Label='cTmp')
    call mma_allocate(sTmp,nTmp,Label='sTmp')
    call mma_allocate(EJ,nConf,Label='EJ')
    call mma_allocate(HEJ,nConf,Label='HEJ')

    call Faroald_Ints(TU,TUVX)
    EJ(:) = Zero
    do J=1,NPCSF
      EJ(IPCSF(J)) = One
      call Faroald_Sigma(nConf,EJ,HEJ,cTmp,sTmp,nTmp,nDetA,nDetB)
      EJ(IPCSF(J)) = Zero
      do I=J,NPCSF
        H(nTri_Elem(I-1)+J) = HEJ(IPCSF(I))
      end do
      H(nTri_Elem(J)) = H(nTri_Elem(J))+ECORE
    end do

    call mma_deallocate(HEJ)
    call mma_deallocate(EJ)
    call mma_deallocate(sTmp)
    call mma_deallocate(cTmp)

  end subroutine IHJ_Faroald_Col

  ! Reference block from the columns H e_J of Lucia or SGUGA (Ref_Sigma) in
  ! the basis of the CI back end, J over the selected CSFs.
  subroutine IHJ_Ref_Col(LabRef,H)
    character(len=*), intent(in) :: LabRef
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J
    real(kind=wp), allocatable :: EJ(:), HEJ(:)

    call mma_allocate(EJ,nConf,Label='EJ')
    call mma_allocate(HEJ,nConf,Label='HEJ')
    EJ(:) = Zero
    do J=1,NPCSF
      EJ(IPCSF(J)) = One
      call Ref_Sigma(LabRef,1,IREFSM,nConf,EJ,HEJ,nTU,TU,nTUVX,TUVX)
      EJ(IPCSF(J)) = Zero
      do I=J,NPCSF
        H(nTri_Elem(I-1)+J) = HEJ(IPCSF(I))
      end do
      H(nTri_Elem(J)) = H(nTri_Elem(J))+ECORE
    end do
    call mma_deallocate(HEJ)
    call mma_deallocate(EJ)

  end subroutine IHJ_Ref_Col

  ! Compare PHP (back end Lab) with the reference block HRef (back end LabRef).
  subroutine Compare_IHJ(Lab,LabRef,HRef)
    character(len=*), intent(in) :: Lab, LabRef
    real(kind=wp), intent(in) :: HRef(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, IJ
    real(kind=wp) :: ErrMax

    ErrMax = Zero
    do I=1,NPCSF
      do J=1,I
        IJ = nTri_Elem(I-1)+J
        ErrMax = max(ErrMax,abs(PHP(IJ)-HRef(IJ)))
        if (abs(PHP(IJ)-HRef(IJ)) > Thr) then
          write(u6,*) 'Mk_IHJ: <I|H|J> of ',Lab,' differs from ',LabRef
          write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IPCSF(I),IPCSF(J)
          write(u6,*) Lab,' =',PHP(IJ)
          write(u6,*) LabRef,' =',HRef(IJ)
          call Abend()
        end if
      end do
    end do
    write(u6,'(1X,5A,I0,A,ES10.3)') 'Mk_IHJ: ',Lab,' verified against ',LabRef,', ',NPCSF,' CSFs, max. error ',ErrMax

  end subroutine Compare_IHJ

  ! Compare PHP with columns H e_J computed by Mk_H_Psi (the sigma vector of
  ! the CI back end, with SIGMA the one of SGUGA) for at most mxCol CSFs J
  ! of the selection, evenly spaced.
  subroutine Verify_IHJ()
    integer(kind=iwp), parameter :: mxCol = 32
    integer(kind=iwp) :: I, J, JJ, nCol, nTmp
    real(kind=wp), allocatable :: EJ(:), HEJ(:), cTmp(:), sTmp(:)
    real(kind=wp) :: Ref, Err, ErrMax
    logical(kind=iwp) :: Sigma_Save

    if (DoFaro) then
      nTmp = Faro_nTmp()
    else
      nTmp = nDet
    end if
    call mma_allocate(cTmp,nTmp,Label='cTmp')
    call mma_allocate(sTmp,nTmp,Label='sTmp')
    call mma_allocate(EJ,nConf,Label='EJ')
    ! CSDTVC may write beyond nConf (cf. David5)
    call mma_allocate(HEJ,max(nConf,nDet),Label='HEJ')

    Sigma_Save = Sigma_on_disk
    nCol = min(NPCSF,mxCol)
    ErrMax = Zero
    do JJ=1,nCol
      J = 1+((JJ-1)*(NPCSF-1))/max(nCol-1,1)
      EJ(:) = Zero
      EJ(IPCSF(J)) = One
      call Mk_H_Psi(iState,IREFSM,nConf,EJ,HEJ,cTmp,sTmp,nTmp,nTU,TU,nTUVX,TUVX)
      do I=1,NPCSF
        Ref = HEJ(IPCSF(I))
        if (I == J) Ref = Ref+ECORE
        Err = abs(PHP(nTri_Elem(max(I,J)-1)+min(I,J))-Ref)
        ErrMax = max(ErrMax,Err)
        if (Err > Thr) then
          write(u6,*) 'Mk_IHJ: <I|H|J> differs from (H e_J)(I)'
          write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IPCSF(I),IPCSF(J)
          write(u6,*) 'Mk_IHJ   =',PHP(nTri_Elem(max(I,J)-1)+min(I,J))
          write(u6,*) 'Mk_H_Psi =',Ref
          call Abend()
        end if
      end do
    end do
    Sigma_on_disk = Sigma_Save
    write(u6,'(1X,A,I0,A,I0,A,ES10.3)') 'Mk_IHJ verified against Mk_H_Psi: ',nCol,' of ',NPCSF,' columns, max. error ',ErrMax

    call mma_deallocate(HEJ)
    call mma_deallocate(EJ)
    call mma_deallocate(sTmp)
    call mma_deallocate(cTmp)

  end subroutine Verify_IHJ

end subroutine Mk_IHJ

!> Excitation level between two configurations given by their doubly (D)
!> and singly (S) occupied orbitals as bit strings: the number of electrons
!> to move, Sum_t max(0,n_t(L)-n_t(R)). H couples only configurations with
!> an excitation level of at most 2.
function Cnf_Exc(DL,SL,DR,SR) result(nExc)
integer(kind=iwp) :: nExc
integer(kind=iwp), intent(in) :: DL, SL, DR, SR

! n_t(L)-n_t(R) = 2 (doubly -> empty), 1 (doubly -> singly, singly -> empty)
nExc = 2*popcnt(iand(DL,not(ior(DR,SR))))+popcnt(iand(DL,SR))+popcnt(iand(SL,not(ior(DR,SR))))

end function Cnf_Exc

!> <aL bL|H|aR bR> between Faroald determinants (alpha and beta strings as
!> bit strings, orbital t in bit t-1, alpha string before the beta string)
!> by the Slater-Condon rules, with the integrals htu and gtuvx (chemists'
!> notation) of Faroald, H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx
!> (E_tu E_vx - delta_uv E_tx), as in sigma_update.
function Faroald_HDet(aL,bL,aR,bR) result(E)
use Constants, only: Half

real(kind=wp) :: E
integer(kind=iwp), intent(in) :: aL, bL, aR, bR
integer(kind=iwp) :: dA, dB, nA, nB, i, j, p, q, r, s, aI

dA = ieor(aL,aR)
dB = ieor(bL,bR)
nA = popcnt(dA)/2
nB = popcnt(dB)/2
E = Zero
if (nA+nB > 2) return

if (nA+nB == 0) then
  do i=1,my_norb
    if (btest(aL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(aL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
        if (btest(bL,j-1)) E = E+gtuvx(i,i,j,j)
      end do
    end if
    if (btest(bL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(bL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
      end do
    end if
  end do
else if (nA == 1 .and. nB == 0) then
  call Single(aL,aR,bR)
else if (nA == 0 .and. nB == 1) then
  call Single(bL,bR,aR)
else if (nA == 2) then
  call Double(aL,aR)
else if (nB == 2) then
  call Double(bL,bR)
else
  ! one alpha (p <- q) and one beta (r <- s) excitation
  p = Orb(iand(dA,aL),1)
  q = Orb(iand(dA,aR),1)
  r = Orb(iand(dB,bL),1)
  s = Orb(iand(dB,bR),1)
  E = Par(aR,p,q)*Par(bR,r,s)*gtuvx(p,q,r,s)
end if

contains

  ! single excitation p <- q in the string xL/xR of one spin, yR the string
  ! of the other spin
  subroutine Single(xL,xR,yR)
    integer(kind=iwp), intent(in) :: xL, xR, yR
    integer(kind=iwp) :: dX, iC, j

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    q = Orb(iand(dX,xR),1)
    iC = iand(xL,xR)
    E = htu(p,q)
    do j=1,my_norb
      if (btest(iC,j-1)) E = E+gtuvx(p,q,j,j)-gtuvx(p,j,j,q)
      if (btest(yR,j-1)) E = E+gtuvx(p,q,j,j)
    end do
    E = Par(xR,p,q)*E

  end subroutine Single

  ! double excitation p,r <- q,s in the string xL/xR of one spin
  subroutine Double(xL,xR)
    integer(kind=iwp), intent(in) :: xL, xR
    integer(kind=iwp) :: dX

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    r = Orb(iand(dX,xL),2)
    q = Orb(iand(dX,xR),1)
    s = Orb(iand(dX,xR),2)
    aI = ibset(ibclr(xR,q-1),p-1)
    E = Par(xR,p,q)*Par(aI,r,s)*(gtuvx(p,q,r,s)-gtuvx(p,s,r,q))

  end subroutine Double

  ! the n-th orbital (bit) set in x
  function Orb(x,n) result(k)
    integer(kind=iwp) :: k
    integer(kind=iwp), intent(in) :: x, n
    integer(kind=iwp) :: m

    m = 0
    do k=1,my_norb
      if (btest(x,k-1)) then
        m = m+1
        if (m == n) return
      end if
    end do
    k = 0

  end function Orb

  ! sign of a+_t a_u acting on the string x: the parity of the number of
  ! orbitals occupied in x strictly between t and u
  function Par(x,t,u) result(iSgn)
    integer(kind=iwp) :: iSgn
    integer(kind=iwp), intent(in) :: x, t, u
    integer(kind=iwp) :: lo, hi

    lo = min(t,u)
    hi = max(t,u)
    iSgn = 1-2*poppar(ibits(x,lo,hi-lo-1))

  end function Par

end function Faroald_HDet

subroutine Mk_Tra_Psi(iState,TU,nTU,CI_Vector,nCSF)

use sguga, only: sg_traci
use general_data, only: nSym

integer(kind=iwp), intent(in) ::  iState,nTU, nCSF
real(kind=wp), intent(in) :: TU(nTU)
real(kind=wp), intent(inout) :: CI_Vector(nCSF)

real(kind=wp), allocatable :: Faroald_Psi(:,:), ctemp(:), TFull(:)
integer(kind=iwp) :: iSym, iOffT, iOffO, kk, pp
integer(kind=iwp), allocatable :: O2L(:)

integer(kind=iwp) :: iCSF
real(kind=wp), allocatable :: CI_Faroald(:)
real(kind=wp), allocatable :: CI_SGUGA(:)
logical(kind=iwp) :: Verify_SG

! SGCI: the vector is in the SGUGA basis, transformed by SGUGA; under SIGMA
! compared with Lucia (vector reordered to the Lucia basis and back).
if (SG_CI()) then
  if (CI_Verify('SIGMA')) then
    call mma_allocate(CI_Faroald,nCSF,Label='CI_Lucia')
    call SG_REORD(iState,CI_St(iState)%STSym,1,nCSF,CI_Vector,CI_Faroald)
    call Lucia_util('TRACI',CI_Vector=CI_Faroald(:),TU=TU(:),nTU=nTU)
  end if
  call SG_TRACI(iState,CI_St(iState)%STSym,TU,nTU,CI_Vector,nCSF)
  if (CI_Verify('SIGMA')) then
    call mma_allocate(CI_SGUGA,nCSF,Label='CI_SGUGA')
    call SG_REORD(iState,CI_St(iState)%STSym,0,nCSF,CI_Faroald,CI_SGUGA)
    if (.not. TraPsi_Verified) then
      write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_Tra_Psi (SGCI): SGUGA CI transformation verified against Lucia, ',nCSF, &
                                     ' CSFs, max. error ',maxval(abs(CI_Vector(:)-CI_SGUGA(:)))
      TraPsi_Verified = .true.
    end if
    if (maxval(abs(CI_Vector(:)-CI_SGUGA(:))) > 1.0e-10_wp) then
      write(u6,*) 'Mk_Tra_Psi (SGCI): SGUGA and Lucia CI transformations differ'
      call Abend()
    end if
    call mma_deallocate(CI_SGUGA)
    call mma_deallocate(CI_Faroald)
  end if
  return
end if

! CI development (MOLCAS_CI_VERIFY=SIGMA): also transform with SGUGA and compare.
Verify_SG = CI_Verify('SIGMA') .and. SG_Avail()
If (Verify_SG) Then
  Call mma_allocate(CI_SGUGA,nCSF,Label='CI_SGUGA')
  CI_SGUGA(:) = CI_Vector(:)
End If

If (DoFaro .and. CI_Verify('SIGMA')) Then
  ! the reference of Lucia, in the Lucia basis
  Call mma_allocate(CI_Faroald,nCSF,Label='CI_Lucia')
  call CI_Basis_To_Lucia(iState,CI_St(iState)%STSym,nCSF,CI_Vector,CI_Faroald)
  Call Lucia_util('TRACI',CI_Vector=CI_Faroald(:),TU=TU(:),nTU=nTU)
End If

If (DoFaro) Then
  ! Faroald in its own basis (citrans), through the determinants
  ! the transformation matrix in symmetry blocks -> NAC x NAC over the
  ! orbitals of Faroald (the levels, Faro_L2O)
  Call mma_allocate(TFull,NAC*NAC,Label='TFull')
  Call mma_allocate(O2L,NAC,Label='O2L')
  do kk=1,NAC
    O2L(Faro_L2O(kk)) = kk
  end do
  TFull(:) = Zero
  iOffT = 0
  iOffO = 0
  do iSym=1,nSym
    do kk=1,NASH(iSym)
      do pp=1,NASH(iSym)
        TFull(O2L(iOffO+pp)+(O2L(iOffO+kk)-1)*NAC) = TU(iOffT+pp+(kk-1)*NASH(iSym))
      end do
    end do
    iOffT = iOffT+NASH(iSym)**2
    iOffO = iOffO+NASH(iSym)
  end do
  Call mma_deallocate(O2L)
  If (ras_active) Then
    ! RAS: in the RAS string classes (faroald_ras)
    Call mma_allocate(Faroald_Psi,max(1,n_ras),1,Label='Faroald_Psi')
    call CITRANS_CSF2RAS(CI_Vector,Faroald_Psi(:,1))
    Call RAS_Traci(TFull,NAC*NAC,Faroald_Psi(:,1))
    call CITRANS_RAS2CSF(Faroald_Psi(:,1),CI_Vector)
  Else
    Call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Faroald_Psi')
    Faroald_PSI(:,:) = Zero
    call CITRANS_CSF2SD(CI_Vector,Faroald_PSI)
    Call traci_faroald(TFull,NAC*NAC,Faroald_Psi)
    call CITRANS_SD2CSF(Faroald_Psi,CI_Vector)
  End If
  Call mma_deallocate(TFull)
  Call mma_deallocate(Faroald_Psi)
  If (CI_Verify('SIGMA')) Then
    Call mma_allocate(ctemp,nCSF,Label='ctemp')
    call CI_Lucia_To_Basis(iState,CI_St(iState)%STSym,nCSF,CI_Faroald,ctemp)
    call Report_Dev('Mk_Tra_Psi (FARO): Faroald CI transformation vs Lucia',maxval(abs(CI_Vector(:)-ctemp(:))),1.0e-10_wp, &
                    TraPsi_Faro_Verified)
    Call mma_deallocate(ctemp)
    Call mma_deallocate(CI_Faroald)
  End If
Else
  ! Note that the CI vector is transformed to Lucia format inside the
  ! routine! This in difference to the other calls to Lucia_util.
  Call Lucia_util('TRACI', &
                  CI_Vector=CI_Vector(:), &
                  TU=TU(:), nTU=nTU)
End If

If (Verify_SG) Then

!  Call SG_TRACI(iState,STSYM,TU,nTU,CI_SGUGA,nCSF)
   Call mma_allocate(ctemp,nCSF,Label='ctemp')
   call CI_Basis_To_File(iState,CI_St(iState)%STSym,nCSF,CI_SGUGA,ctemp)
   Call SG_TRACI(iState,CI_St(iState)%STSym,TU,nTU,ctemp,nCSF)
   call CI_File_To_Basis(iState,CI_St(iState)%STSym,nCSF,ctemp,CI_SGUGA)
   Call mma_deallocate(ctemp)
   If (.not. TraPsi_Verified) Then
      Write (u6,'(1X,A,I0,A,ES10.3)') 'Mk_Tra_Psi: SGUGA CI transformation verified against the CI back end, ',nCSF, &
                                     ' CSFs, max. error ',maxval(abs(CI_Vector(:)-CI_SGUGA(:)))
      TraPsi_Verified = .true.
   End If

   Do iCSF = 1, nCSF
      If (Abs(CI_Vector(iCSF)-CI_SGUGA(iCSF))/nCSF>1.0E-12_wp) Then
         Write (6,*) 'Mk_Tra_Psi(SGUGA) error, iCSF=',iCSF
         Write (6,*) 'CI_Vector(iCSF) =',CI_Vector(iCSF)
         Write (6,*) 'CI_SGUGA(iCSF)=',CI_SGUGA(iCSF)
         Call Abend()
      End If
   End Do
   Call mma_deallocate(CI_SGUGA)
End If

End subroutine Mk_Tra_Psi

End module CI_Interfaces
