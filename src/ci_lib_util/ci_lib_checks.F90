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

!> The verification layer of the CI interface (MOLCAS_CI_VERIFY, see
!> CI_Verify_Options) and its helpers: which library runs (CI_SG_CI, CI_Faro_CI,
!> CI_SG_Avail), conversions of CI vectors to the Lucia and file orders for the
!> references, the reference sigma vectors (CI_Chk_Ref_Sigma), the SGUGA wrappers
!> for GAS (relaxed space), the flags of the once-per-run checks and the
!> timers of the CI operations. It uses only the core (CI_Lib_Core).
module CI_Lib_Checks

use CI_SG_Lucia_Bridge, only: CI_SG_ReOrd
use Lucia_Interface, only: Lucia_Util
use Faroald, only: Faro_sort
use CI_Lib_Defs, only: CI_Lib_Selected, CI_Lib_Faroald, CI_Lib_Is, CI_Lib_SGUGA
use CI_Lib_Core, only: CI_Lucia_Active, CI_HEXS_Keep_SG, CI_HEXS_Proj, CI_SG_Rel_Embed, CI_SG_Rel_Project
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use definitions, only: wp, iwp, u6
use sguga, only: CIS, SGS, EXS
use SGUGA_Sigma, only: SG_Work_Need, sg_h_psi, sg_h_psi_ras, sg_hdiag_fast, sg_hdiag_ref
use SGUGA_Densities, only: sg_spin_pdm, sg_trans_pdm, sg_two_pdm
use CI_Input, only: iDoGAS
use CI_Lib_Info, only: CI_Lib_Space

implicit none
private

public :: CI_Basis, CI_Basis_To_File, CI_Basis_To_Lucia, CI_Faro_nDet_RAS, CI_File_Index_To_Basis, CI_File_To_Basis, &
          CI_Lucia_To_Basis, CI_Timer, CI_Work_Need, CI_Faro_CI, CI_Chk_HPsi_Faro, CI_Chk_HPsi_FaroRAS, &
          CI_Chk_HPsi_RAS, CI_Chk_HPsi_SG, CI_Chk_IHI_Avg, CI_Chk_IHI_Faro_SG, CI_Chk_IHI_SG, CI_Chk_IHJ, CI_Tim_nOp, &
          CI_Chk_PDM_Faro, CI_Chk_PDM_SG, CI_Chk_Ref_Sigma, CI_Chk_SD_ASTA, CI_SG_Avail, CI_SG_CI, CI_SG_HDiag_G, CI_SG_HPsi, &
          CI_SG_Spin_Density_ASTA, CI_SG_Spin_PDM_G, CI_SG_Trans_PDM_G, CI_SG_Two_PDM_G, CI_Chk_T1DM_Faro, CI_Chk_T1DM_SG, &
            CI_Tim_Calls, &
          CI_Tim_CPU, CI_Tim_HPsi, CI_Tim_IHI, CI_Tim_IHJ, CI_Tim_Label, CI_Tim_Off, CI_Tim_On, CI_Tim_PDM, CI_Tim_Print, &
            CI_Tim_T1DM, CI_Tim_TraPsi, CI_Tim_Wall, &
          CI_Chk_TraPsi_Faro, CI_Chk_TraPsi

! CI_Lib_Block is verified against CI_Lib_Sigma only once per run.
logical(kind=iwp) :: CI_Chk_IHJ = .false.
! sg_h_psi_ras is compared with sg_h_psi (CAS) only once per run; for RAS
! the difference to the CI back end is printed once (checked every time).
logical(kind=iwp) :: CI_Chk_HPsi_RAS = .false.
! The SGUGA CI transformation (CI_Lib_Transform, SIGMA) is checked on every call,
! its maximum error printed once.
logical(kind=iwp) :: CI_Chk_TraPsi = .false.
! SGCI: the SGUGA sigma vector compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: CI_Chk_HPsi_SG = .false.
! SGCI: the SGUGA densities compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: CI_Chk_PDM_SG = .false.
! SGCI: the SGUGA diagonal checked against sg_hdiag_ref (SIGMA), printed once.
logical(kind=iwp) :: CI_Chk_IHI_SG = .false.
! SGCI: the SGUGA configuration-average diagonal checked against Lucia, printed once.
logical(kind=iwp) :: CI_Chk_IHI_Avg = .false.
! SGCI: the SGUGA transition densities compared with Lucia (SIGMA), printed once.
logical(kind=iwp) :: CI_Chk_T1DM_SG = .false.
! FARO: the Faroald steps verified against Lucia and SGUGA (SIGMA), printed once
logical(kind=iwp) :: CI_Chk_HPsi_Faro = .false., CI_Chk_PDM_Faro = .false., CI_Chk_T1DM_Faro = .false., &
                     CI_Chk_TraPsi_Faro = .false., CI_Chk_IHI_Faro_SG = .false.
! FARO with RAS: the string-driven sigma (Faroald_RAS) checked once against
! the sigma of the full determinant space (SIGMA)
logical(kind=iwp) :: CI_Chk_HPsi_FaroRAS = .false.
! SIGMA: the spin density of the ASTA triplet operator checked, printed once
logical(kind=iwp) :: CI_Chk_SD_ASTA = .false.

! CPU and wall time and number of calls of the operations of the CI
! interface, for every back end (CI_Timer: 'Init' resets, 'Print' prints)
integer(kind=iwp), parameter :: CI_Tim_HPsi = 1, CI_Tim_IHI = 2, CI_Tim_IHJ = 3, CI_Tim_PDM = 4, CI_Tim_T1DM = 5, &
  CI_Tim_TraPsi = 6, CI_Tim_nOp = 6
character(len=*), parameter :: CI_Tim_Label(CI_Tim_nOp) = ['H|Psi> (CI_Lib_Sigma)          ','diagonal (CI_Lib_Diag)         ', &
                                                    'explicit block (CI_Lib_Block)  ','densities (CI_Lib_Dens)        ', &
                                                    'transition (CI_Lib_Dens)       ','CI transform (CI_Lib_Transform)']
real(kind=wp) :: CI_Tim_CPU(CI_Tim_nOp) = Zero, CI_Tim_Wall(CI_Tim_nOp) = Zero
integer(kind=iwp) :: CI_Tim_Calls(CI_Tim_nOp) = 0

contains

!> SGCI: the whole CI (diagonal, explicit subspace, sigma vectors, density
!> matrices) is done by SGUGA in its own split-graph CSF basis, without
!> Lucia (RASSCF keyword SGCI: SGUGA is the selected CI library, CI_Lib_Defs).
function CI_SG_CI() result(SGMode)
logical(kind=iwp) :: SGMode
SGMode = CI_Lib_Is(CI_Lib_SGUGA)
end function CI_SG_CI

!> FARO: the whole CI is done by Faroald in the configuration order of
!> Faroald_CITrans (RASSCF keyword FARO: Faroald is the selected CI library).
function CI_Faro_CI() result(FaroMode)
logical(kind=iwp) :: FaroMode
FaroMode = CI_Lib_Is(CI_Lib_Faroald)
end function CI_Faro_CI

!> FARO: the number of determinants of the RAS/GAS space of Faroald (all
!> symmetries; nela/nelb electrons, o(k) orbitals in partition k of np,
!> mn(k) to mx(k) electrons in the partitions 1..k), before any setup, for
!> the memory estimate of RASSCF.
function CI_Faro_nDet_RAS(nela,nelb,np,o,mn,mx) result(n)
use Faroald, only: Faro_RAS_count_dets
integer(kind=iwp), intent(in) :: nela, nelb, np, o(np), mn(np), mx(np)
real(kind=wp) :: n
n = Faro_RAS_count_dets(nela,nelb,np,o,mn,mx)
end function CI_Faro_nDet_RAS

!> The memory (words) the active CI back ends need besides the CI vectors of
!> the caller, for the Davidson vector store (Ini_David) to leave free; the
!> programs get it only from here, never from the CI libraries. Lucia (if
!> initialized, also as the reference of MOLCAS_CI_VERIFY=SIGMA): its own
!> estimate; SGUGA (SGCI): the blocks of its sigma and 2-RDM (for GAS in the
!> relaxed space, plus the two vectors there); Faroald sizes its work arrays
!> from the free memory and needs nothing reserved.
function CI_Work_Need() result(nNeed)
use lucia_data, only: Memory_Needed_Lucia
integer(kind=iwp) :: nNeed
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: iR, LSM
nNeed = 0
if (CI_Lucia_Active()) nNeed = nNeed+Memory_Needed_Lucia
if (.not. CI_SG_CI()) return
LSM = CI_Lib_Space(iState)%STSym
iR = CI_Lib_Space(iState)%iRel
if (iR > 0) then
  nNeed = nNeed+SG_Work_Need(SGS(iR),CIS(iR),LSM,.true.)+2*CIS(iR)%NCSF(LSM)
else
  nNeed = nNeed+SG_Work_Need(SGS(iState),CIS(iState),LSM,CI_Lib_Space(iState)%nRsPrt /= 1)
end if
end function CI_Work_Need


!> The CSF basis of the CI vectors of RASSCF = the selected CI library
!> (CI_Lib_Defs): 1 = Lucia, 2 = the split-graph basis of SGUGA (SGCI; the
!> order of the files), 3 = the configuration order of Faroald (Faroald_CITrans;
!> keyword FARO).
function CI_Basis() result(iBasis)
integer(kind=iwp) :: iBasis
iBasis = CI_Lib_Selected
end function CI_Basis

!> A CI vector from the basis of the CI back end to the split-graph order of
!> the JOBIPH/HDF5 files (and of SGUGA).
subroutine CI_Basis_To_File(iState,LSM,n,Cin,Cout)
integer(kind=iwp), intent(in) :: iState, LSM, n
real(kind=wp), intent(in) :: Cin(n)
real(kind=wp), intent(out) :: Cout(n)
select case (CI_Basis())
  case (1)
    call CI_SG_ReOrd(iState,LSM,0,n,Cin,Cout)
  case (2)
    Cout(:) = Cin(:)
  case (3)
    call Faro_sort('O',Cin,Cout)
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
    call CI_SG_ReOrd(iState,LSM,1,n,Cin,Cout)
  case (2)
    Cout(:) = Cin(:)
  case (3)
    call Faro_sort('C',Cin,Cout)
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
  call CI_SG_ReOrd(iState,LSM,1,n,Tmp,Cout)
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
  call CI_SG_ReOrd(iState,LSM,0,n,Cin,Tmp)
  call CI_File_To_Basis(iState,LSM,n,Tmp,Cout)
  call mma_deallocate(Tmp)
end if
end subroutine CI_Lucia_To_Basis

!> Reference sigma vectors (CI development, MOLCAS_CI_VERIFY=SIGMA): Sout =
!> H Cin, both in the basis of the CI back end, computed by Lucia (Lab =
!> 'Lucia', the vector reordered to the Lucia basis and back) or by SGUGA
!> (Lab = 'SGUGA', through the split-graph order of the files).
subroutine CI_Chk_Ref_Sigma(Lab,iState,LSM,n,Cin,Sout,nTU,TU,nTUVX,TUVX)
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
    nD = CI_Lib_Space(iState)%nDet(LSM)
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
    call CI_HEXS_Proj(CI_HEXS_Keep_SG,C1)
    call CI_SG_HPsi(iState,C1,n,LSM,C2,TUVX,nTUVX,TU,nTU,CI_Lib_Space(iState)%nRsPrt /= 1)
    call CI_HEXS_Proj(CI_HEXS_Keep_SG,C2)
    call CI_File_To_Basis(iState,LSM,n,C2,Sout)
  case default
    write(u6,*) 'CI_Chk_Ref_Sigma: unknown back end ',Lab
    call Abend()
end select
call mma_deallocate(C2)
call mma_deallocate(C1)
end subroutine CI_Chk_Ref_Sigma

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


!> SGUGA sigma vector of state iState: sg_h_psi_ras (RAS, or Use_RAS) or
!> sg_h_psi (CAS); for GAS sg_h_psi_ras in the relaxed space.
subroutine CI_SG_HPsi(iState,Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU,Use_RAS)
integer(kind=iwp), intent(in) :: iState, n, LSM, nTUVX, nTU
real(kind=wp), intent(in) :: Psi(n), TUVX(nTUVX), TU(nTU)
real(kind=wp), intent(out) :: Sigma(n)
logical(kind=iwp), intent(in) :: Use_RAS
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: PR(:), SR(:)
iR = CI_Lib_Space(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(PR,nR,Label='PR')
  call mma_allocate(SR,nR,Label='SR')
  call CI_SG_Rel_Embed(iState,LSM,n,Psi,PR)
  call sg_h_psi_ras(SGS(iR),CIS(iR),EXS(iR),PR,nR,LSM,SR,TUVX,nTUVX,TU,nTU)
  call CI_SG_Rel_Project(iState,LSM,n,SR,Sigma)
  call mma_deallocate(SR)
  call mma_deallocate(PR)
else if (Use_RAS) then
  call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU)
else
  call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),Psi,n,LSM,Sigma,TUVX,nTUVX,TU,nTU)
end if
end subroutine CI_SG_HPsi

!> SGUGA 2-RDM (sg_two_pdm) of state iState; for GAS in the relaxed space.
subroutine CI_SG_Two_PDM_G(iState,C,n,LSM,P,PA,nP)
integer(kind=iwp), intent(in) :: iState, n, LSM, nP
real(kind=wp), intent(in) :: C(n)
real(kind=wp), intent(out) :: P(nP), PA(nP)
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: CR(:)
iR = CI_Lib_Space(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(CR,nR,Label='CR')
  call CI_SG_Rel_Embed(iState,LSM,n,C,CR)
  call sg_two_pdm(SGS(iR),CIS(iR),EXS(iR),CR,nR,LSM,P,PA,nP)
  call mma_deallocate(CR)
else
  call sg_two_pdm(SGS(iState),CIS(iState),EXS(iState),C,n,LSM,P,PA,nP)
end if
end subroutine CI_SG_Two_PDM_G

!> SGUGA spin density (sg_spin_pdm, D the 1-RDM) of state iState; for GAS
!> in the relaxed space.
subroutine CI_SG_Spin_PDM_G(iState,C,n,LSM,D,SD,nD)
integer(kind=iwp), intent(in) :: iState, n, LSM, nD
real(kind=wp), intent(in) :: C(n), D(nD)
real(kind=wp), intent(out) :: SD(nD)
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: CR(:)
iR = CI_Lib_Space(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(CR,nR,Label='CR')
  call CI_SG_Rel_Embed(iState,LSM,n,C,CR)
  call sg_spin_pdm(SGS(iR),CIS(iR),EXS(iR),CR,nR,LSM,D,SD,nD)
  call mma_deallocate(CR)
else
  call sg_spin_pdm(SGS(iState),CIS(iState),EXS(iState),C,n,LSM,D,SD,nD)
end if
end subroutine CI_SG_Spin_PDM_G

!> SGUGA transition 1-RDM and spin 1-TDM (sg_trans_pdm) of state iState;
!> for GAS in the relaxed space.
subroutine CI_SG_Trans_PDM_G(iState,B,K,n,LSM,nA,TS,TSS,Spin)
integer(kind=iwp), intent(in) :: iState, n, LSM, nA
real(kind=wp), intent(in) :: B(n), K(n)
real(kind=wp), intent(out) :: TS(nA*nA), TSS(nA*nA)
logical(kind=iwp), intent(in) :: Spin
integer(kind=iwp) :: iR, nR
real(kind=wp), allocatable :: BR(:), KR(:)
iR = CI_Lib_Space(iState)%iRel
if (iR > 0) then
  nR = CIS(iR)%NCSF(LSM)
  call mma_allocate(BR,nR,Label='BR')
  call mma_allocate(KR,nR,Label='KR')
  call CI_SG_Rel_Embed(iState,LSM,n,B,BR)
  call CI_SG_Rel_Embed(iState,LSM,n,K,KR)
  call sg_trans_pdm(SGS(iR),CIS(iR),EXS(iR),BR,KR,nR,LSM,nA,TS,TSS,Spin)
  call mma_deallocate(KR)
  call mma_deallocate(BR)
else
  call sg_trans_pdm(SGS(iState),CIS(iState),EXS(iState),B,K,n,LSM,nA,TS,TSS,Spin)
end if
end subroutine CI_SG_Trans_PDM_G

!> Exact SGUGA diagonal <I|H|I> of state iState (Mode 'FAST': sg_hdiag_fast,
!> 'REF': sg_hdiag_ref); for GAS in the relaxed space.
subroutine CI_SG_HDiag_G(Mode,iState,n,LSM,TUVX,nTUVX,TU,nTU,Diag)
character(len=*), intent(in) :: Mode
integer(kind=iwp), intent(in) :: iState, n, LSM, nTUVX, nTU
real(kind=wp), intent(in) :: TUVX(nTUVX), TU(nTU)
real(kind=wp), intent(out) :: Diag(n)
integer(kind=iwp) :: iR, iS, nR
real(kind=wp), allocatable :: DR(:)
iR = CI_Lib_Space(iState)%iRel
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
  call CI_SG_Rel_Project(iState,LSM,n,DR,Diag)
else
  Diag(:) = DR(:)
end if
call mma_deallocate(DR)
end subroutine CI_SG_HDiag_G

!> SGUGA covers the CI space of RASSCF: CAS and RAS directly, GAS through
!> the relaxed space (set up by CI_Lib_Initialize; HEXS/DEXS are masks).
function CI_SG_Avail() result(Avail)
integer(kind=iwp), parameter :: iState = 1
logical(kind=iwp) :: Avail
Avail = (.not. iDoGAS) .or. (CI_Lib_Space(iState)%iRel > 0)
end function CI_SG_Avail

!> One-particle spin density SD(iTri(p,q)) = <Psi|a+_pa a_qa - a+_pb a_qb|Psi>
!> (levels, triangular packing as sg_spin_pdm) of the CI vector CI in the
!> split-graph order, from the ASTA triplet operator (SG_Tpq_Psi): the
!> segment values of the triplet operator, no two-particle intermediates.
subroutine CI_SG_Spin_Density_ASTA(iState,LSM,n,CI,nTri,SD)
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
end subroutine CI_SG_Spin_Density_ASTA

!> Start the timer of an operation (T: CPU and wall time at the start).
subroutine CI_Tim_On(T)
real(kind=wp), intent(out) :: T(2)
call CWTime(T(1),T(2))
end subroutine CI_Tim_On

!> Stop the timer of operation k started at T.
subroutine CI_Tim_Off(k,T)
integer(kind=iwp), intent(in) :: k
real(kind=wp), intent(in) :: T(2)
real(kind=wp) :: C, W
call CWTime(C,W)
CI_Tim_CPU(k) = CI_Tim_CPU(k)+C-T(1)
CI_Tim_Wall(k) = CI_Tim_Wall(k)+W-T(2)
CI_Tim_Calls(k) = CI_Tim_Calls(k)+1
end subroutine CI_Tim_Off

subroutine CI_Timer(Opt)
use lucia_data, only: TDENSI, TSIGMA

character(Len=*) Opt

Select case (Opt(1:3))
Case ('Ini')
   !SVC: lucia timers
   tsigma(:) = Zero
   tdensi(:) = Zero
   CI_Tim_CPU(:) = Zero
   CI_Tim_Wall(:) = Zero
   CI_Tim_Calls(:) = 0
Case ('Pri')
  call CI_Tim_Print()
  If (CI_Faro_CI()) Then
  Else If (CI_Lucia_Active()) Then
    call PrtTim()
  End If
Case Default
  Write (u6,*) 'CI_Timer called with illegal option:',Opt
  Call Abend()
End Select

end subroutine CI_Timer

!> The timings of the CI interface (CPU, wall, calls per operation).
subroutine CI_Tim_Print()
integer(kind=iwp) :: k
character(len=8) :: BackEnd
if (CI_SG_CI()) then
  BackEnd = 'SGUGA'
else if (CI_Faro_CI()) then
  BackEnd = 'Faroald'
else
  BackEnd = 'Lucia'
end if
write(u6,*)
write(u6,'(1X,A,A)') 'CI interface timings, back end ',trim(BackEnd)
write(u6,'(1X,A32,3X,A12,A12,A10)') 'operation','CPU (s)','wall (s)','calls'
do k=1,CI_Tim_nOp
  if (CI_Tim_Calls(k) > 0) write(u6,'(1X,A32,3X,2F12.3,I10)') CI_Tim_Label(k),CI_Tim_CPU(k),CI_Tim_Wall(k),CI_Tim_Calls(k)
end do
write(u6,'(1X,A32,3X,2F12.3)') 'total',sum(CI_Tim_CPU),sum(CI_Tim_Wall)
end subroutine CI_Tim_Print

end module CI_Lib_Checks
