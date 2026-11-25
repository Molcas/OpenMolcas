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
use CI_Lib_Checks
use CI_Library, only: Checked_Sigma => CI_Lib_Sigma
use SG_Lucia_Bridge, only: SG_ReOrd, SG_Verify_Lucia_Walk
use Lucia_Interface, only: Lucia_Util
use faroald, only: faroald_full, my_norb, sigma_update, htu, gtuvx, ndeta, ndetb
use citrans, only: citrans_csf2ras, citrans_csf2sd, citrans_ras2csf, citrans_sd2csf, citrans_sort
use faroald_ras, only: n_ras, ras_active, ras_sigma
! FARO: Faroald works with the active orbitals in the order of the levels
! of the split graph (by partition, within each by symmetry), so that its
! CSFs couple the spins in the order of the file basis. Faro_L2O(lev) is
! the active orbital (symmetry-blocked order of TU/TUVX) of level lev, set
! up by the one Faroald setup of CI_Operations; the densities are in level
! order (as Lucia's and SGUGA's).
use CI_Operations, only: CI_Lucia_Active, HEXS_Keep_SG, HEXS_Proj
use CI_Operations, only: CI_Lib_Block, CI_Lib_Cnf_Layout, CI_Lib_nCnf, CI_Lib_SGEX, Faro_L2O => Faro_Lev, Lucia1_Block, &
                         Lucia1_Block_SG, Rel_Map, Report_Dev, SG_Rel_Map
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use definitions, only: wp, iwp, u6
use sguga, only: CIS, SGS, SGCnfList, SG_Explicit_Block
use CI_Verify_Options, only: CI_Verify
use CI_Info, only: CI_St

#include "intent.fh"

implicit none
Private

Public :: Mk_H_Psi_Det, Mk_IHJ
Public :: CI_Basis, CI_Basis_To_File, CI_Faro_nDet_RAS, CI_File_Index_To_Basis, CI_File_To_Basis, &
          CI_Lucia_Active, CI_Timer, CI_Verify, CI_Work_Need, Faro_CI, Faro_nTmp, SG_CI

! CI_Verify (module CI_Verify_Options) is re-exported, see there for the
! run-time verification options (MOLCAS_CI_VERIFY).

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************


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

if (Faro_CI()) then

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


!***********************************************************************************************************************************
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!***********************************************************************************************************************************





!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Explicit CI Hamiltonian <I|H|J> over the CSFs of the configurations
!> IPCNF(1:NPCNF) of symmetry IREFSM (configuration indices as in
!> CI_Lib_Cnf_Layout), packed as a lower triangle: the CSFs are ordered
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
!>  - Faroald (FARO): for each pair of configurations the determinant
!>    block (Slater-Condon rules, Faroald_HDet) transformed with the spin
!>    tables of citrans, C_L^T H C_R.
!>  - SGUGA (RASSCF keyword SGEX; no GAS, single RAS partition): <I|H|J>
!>    for each pair of CSFs from the occupation difference (SG_MatEl_H),
!>    pairs of configurations more than a double excitation apart skipped.
!> With MOLCAS_CI_VERIFY=SIGMA, once per run, the block in use is compared
!> with the other back ends (the Faroald block also with columns H e_J of
!> the Faroald sigma vector) and checked against Mk_H_Psi; the CPU time of
!> each back end is reported.
subroutine Mk_IHJ(PHP,nPHP,IPCNF,NPCNF,IREFSM,ECORE,TU,nTU,TUVX,nTUVX)
use Index_Functions, only: nTri_Elem
use Lucia_data, only: Sigma_on_disk
use Constants, only: One

integer(kind=iwp), intent(in) :: nPHP, NPCNF, IPCNF(NPCNF), IREFSM, nTU, nTUVX
real(kind=wp), intent(out) :: PHP(nPHP)
real(kind=wp), intent(in) :: ECORE, TU(nTU), TUVX(nTUVX)
integer(kind=iwp) :: NPCSF
integer(kind=iwp), allocatable :: IPCSF(:), nSelCnf(:)
real(kind=wp) :: T0, T1, W
real(kind=wp), allocatable :: PHP_Ref(:)
real(kind=wp), parameter :: Thr = 1.0e-10_wp
! the CI state of RASSCF (SGS, CIS, EXS)
integer(kind=iwp), parameter :: iState = 1
logical(kind=iwp) :: SG_OK, UseSG, Verify
! the CSFs of the CI symmetry
integer(kind=iwp) :: nConf
character(len=8) :: Lab

! CSF indices (in the CSF basis of symmetry IREFSM) of the selection.
real(kind=wp) :: TimT(2)
call Tim_On(TimT)
nConf = CI_St(iState)%nCSF(IREFSM)
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
  call Lib_Block(PHP)
  if (CI_Verify('SIGMA') .and. (.not. IHJ_Verified)) then
    IHJ_Verified = .true.
    call Verify_IHJ()
  end if
  call mma_deallocate(IPCSF)
  call mma_deallocate(nSelCnf)
  call Tim_Off(Tim_IHJ,TimT)
  return
end if

! SGUGA is available without GAS and for a single RAS partition (SGS is not
! set up with GAS, so test iDoGAS first).
SG_OK = SG_Avail()
! SGEX: the block of the core is SGUGA's for RASSCF's Lucia (CI_Lib_SGEX)
UseSG = CI_Lib_SGEX(iState)
Verify = CI_Verify('SIGMA') .and. (.not. IHJ_Verified)
! the walks of the Lucia CSFs, once per run
if (Verify .and. SG_OK .and. (.not. Faro_CI())) call SG_Verify_Lucia_Walk(iState,IREFSM,nConf)

call CWTime(T0,W)
if (UseSG) then
  Lab = 'SGUGA'
  call Lib_Block(PHP)
else if (Faro_CI()) then
  Lab = 'Faroald'
  call Lib_Block(PHP)
else
  Lab = 'Lucia'
  call Lib_Block(PHP)
end if
call CWTime(T1,W)

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare the block, once per run,
! with the other back ends and with exact columns H e_J from Mk_H_Psi.
if (Verify) then
  write(u6,'(1X,3A,I0,A,F10.3,A)') 'Mk_IHJ: ',trim(Lab),' block, ',NPCSF,' CSFs, CPU ',T1-T0,' s'
  call mma_allocate(PHP_Ref,nTri_Elem(NPCSF),Label='PHP_Ref')
  ! (with FARO the block is in the basis of Faroald: the references of
  ! Lucia and SGUGA are their columns H e_J in that basis, Ref_Sigma)
  if (UseSG .and. (.not. Faro_CI())) call Ref_Block('Lucia')
  if (Faro_CI()) then
    call Ref_Block('Faroald H e_J')
    call Ref_Block('Lucia H e_J')
    if (SG_OK) call Ref_Block('SGUGA H e_J')
  end if
  if (SG_OK) then
    if ((.not. UseSG) .and. (.not. Faro_CI())) call Ref_Block('SGUGA')
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

call Tim_Off(Tim_IHJ,TimT)
contains

  ! The block of the selected library (CI_Lib_Block), for the
  ! configurations IPCNF in its configuration order; its CSFs must be
  ! those of the selection (Sel_CSFs)
  subroutine Lib_Block(HB)
    real(kind=wp), intent(out) :: HB(*)
    integer(kind=iwp) :: nPB
    integer(kind=iwp), allocatable :: IdxB(:)
    call mma_allocate(IdxB,max(1,NPCSF),Label='IdxB')
    call CI_Lib_Block(iState,IREFSM,nTU,TU,nTUVX,TUVX,ECORE,NPCNF,IPCNF,NPCSF,nPB,IdxB,nTri_Elem(NPCSF),HB)
    if ((nPB /= NPCSF) .or. any(IdxB(1:NPCSF) /= IPCSF(1:NPCSF))) then
      write(u6,*) 'Mk_IHJ: the CSFs of CI_Lib_Block differ from the selection'
      call Abend()
    end if
    call mma_deallocate(IdxB)
  end subroutine Lib_Block

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
        ! (with SGEX: the Lucia block of the core)
        call Lucia1_Block(IREFSM,nTU,TU,nTUVX,TUVX,ECORE,NPCNF,IPCNF,NPCSF,PHP_Ref)
      case ('Faroald H e_J')
        call IHJ_Faroald_Col(PHP_Ref)
      case ('Lucia H e_J')
        call IHJ_Ref_Col('Lucia',PHP_Ref)
      case ('SGUGA H e_J')
        call IHJ_Ref_Col('SGUGA',PHP_Ref)
      case ('SGUGA')
        call Lucia1_Block_SG(iState,IREFSM,nTU,TU,nTUVX,TUVX,ECORE,NPCNF,nSelCnf,NPCSF,IPCSF,PHP_Ref)
    end select
    call CWTime(T3,W2)
    write(u6,'(1X,3A,F10.3,A)') 'Mk_IHJ: ',LabRef,' block, CPU ',T3-T2,' s'
    call Compare_IHJ(trim(Lab),LabRef,PHP_Ref)

  end subroutine Ref_Block

  ! IPCSF: the CSFs of the selected configurations, in the order of PHP;
  ! nSelCnf: the number of CSFs of each selected configuration.
  subroutine Sel_CSFs()
    integer(kind=iwp) :: iCnf, iOff, kCnf, I, ICSF, nCnfL
    integer(kind=iwp), allocatable :: nCSFCnf(:), iCSF1(:), iCSFCnf(:)

    nCnfL = CI_Lib_nCnf(iState,IREFSM)
    call mma_allocate(nCSFCnf,nCnfL,Label='nCSFCnf')
    call mma_allocate(iCSFCnf,nConf,Label='iCSFCnf')
    call CI_Lib_Cnf_Layout(iState,IREFSM,nCnfL,nCSFCnf,nConf,iCSFCnf)
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
    call mma_allocate(nCSFL,CI_Lib_nCnf(iState,IREFSM),Label='nCSFL')
    call mma_allocate(GtoL,nConf,Label='GtoL')
    call CI_Lib_Cnf_Layout(iState,IREFSM,size(nCSFL),nCSFL,nConf,GtoL)
    call mma_deallocate(GtoL)
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
    integer(kind=iwp) :: I, J, JJ, nCol
    real(kind=wp), allocatable :: EJ(:), HEJ(:)
    real(kind=wp) :: Ref, Err, ErrMax
    logical(kind=iwp) :: Sigma_Save

    call mma_allocate(EJ,nConf,Label='EJ')
    call mma_allocate(HEJ,nConf,Label='HEJ')

    Sigma_Save = Sigma_on_disk
    nCol = min(NPCSF,mxCol)
    ErrMax = Zero
    do JJ=1,nCol
      J = 1+((JJ-1)*(NPCSF-1))/max(nCol-1,1)
      EJ(:) = Zero
      EJ(IPCSF(J)) = One
      call Checked_Sigma(iState,0,IREFSM,IREFSM,.false.,.true.,nTU,TU,nTUVX,TUVX,nConf,EJ,nConf,HEJ)
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

  end subroutine Verify_IHJ

end subroutine Mk_IHJ

End module CI_Interfaces
