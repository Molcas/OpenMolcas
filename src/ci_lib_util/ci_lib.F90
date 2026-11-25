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

!> The CI library of the CI interface, the module the programs use: the
!> generic CI operations CI_Lib_* of the core (CI_Lib_Core, see there),
!> with the checks of MOLCAS_CI_VERIFY (CI_Lib_Checks) around those that
!> are verified. The checks apply to the CI state set up by
!> CI_Lib_Initialize from the CI input (CI_Lib_Input_State); a state the
!> program set up itself gets the core operations unchanged.
module CI_Lib

! (for the programs: the timers, work-space estimate and conversions of the
! checks module, and the verification options)
use CI_Lib_Checks, only: CI_Faro_nDet_RAS, CI_File_Index_To_Basis, CI_Timer, CI_Work_Need
use CI_Verify_Options, only: CI_Verify
use CI_Lib_Core, Raw_Block => CI_Lib_Block, Raw_Dens => CI_Lib_Dens, Raw_Diag => CI_Lib_Diag, Raw_Sigma => CI_Lib_Sigma, &
                   Raw_Sigma_Det => CI_Lib_Sigma_Det, Raw_Transform => CI_Lib_Transform

implicit none
private :: Block_Col_Ref, Block_Columns, Block_Compare, Block_Eigenvalues, Block_Ref, Check_Block, Check_Diag_Faro, &
           Check_Diag_Faro_SG, Check_Diag_Lucia, Check_Diag_SG, Check_PDMs, Check_Sigma, Check_T1DM, CheckSum, CI_Cnf_Rank, &
           Diag_To_Lucia, Exact_Diag_Lucia, Lucia_Avg_Diag, Lucia_Sigma_Ref, Raw_Block, Raw_Dens, Raw_Diag, Raw_Sigma, &
           Raw_Sigma_Det, Raw_Transform, Sel_PHP, Verify_SG_Explicit

public :: CI_Faro_nDet_RAS, CI_File_Index_To_Basis, CI_Lib_Block, CI_Lib_Dens, CI_Lib_Diag, CI_Lib_Sigma, CI_Lib_Sigma_Det, &
          CI_Lib_Transform, CI_Timer, CI_Verify, CI_Work_Need

contains

!> The 1- and 2-particle (transition) densities of the core (CI_Lib_Dens of
!> CI_Lib_Core, see there for the layouts). Under MOLCAS_CI_VERIFY=SIGMA a
!> one-body transition density of one symmetry (spin-free, not packed) is
!> compared, together with its spin part, with Lucia and SGUGA (those not
!> in use), and so are the packed densities of a state (D, SD, P, PA), the
!> spin density also with the ASTA triplet operator.
subroutine CI_Lib_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Packed)

  use CI_Lib_Checks, only: CI_Tim_Off, CI_Tim_On, CI_Tim_PDM, CI_Tim_T1DM
  use CI_Verify_Options, only: CI_Verify
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, I12, iOpSpin, iLSym, iRSym, nL, nR, n1, n2
  real(kind=wp), intent(in) :: L(nL), R(nR)
  real(kind=wp), intent(out) :: Rho1(n1), Rho2(n2)
  logical(kind=iwp), intent(in), optional :: Packed
  logical(kind=iwp) :: Pk
  real(kind=wp) :: TimT(2)

  Pk = .false.
  if (present(Packed)) Pk = Packed
  call CI_Tim_On(TimT)
  call Raw_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Packed)
  if (CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (.not. Pk) .and. (I12 == 1) .and. (iOpSpin == 0) .and. &
      (iLSym == iRSym) .and. (nL == nR)) call Check_T1DM(iState,iLSym,nL,L,R,n1,Rho1,n2)
  if (CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. Pk .and. (iOpSpin == 0)) &
    call Check_PDMs(iState,iLSym,nL,L,I12,n1,Rho1,n2,Rho2)
  if (Pk) then
    call CI_Tim_Off(CI_Tim_PDM,TimT)
  else
    call CI_Tim_Off(CI_Tim_T1DM,TimT)
  end if

end subroutine CI_Lib_Dens

!> CI development (SIGMA): the packed densities of the state C (D in Rho1;
!> P, PA in Rho2 for I12 = 2; SD computed here) against the libraries not
!> in use: Lucia (vector in its basis; its spin density only for S > 0) and
!> SGUGA (split-graph order; not for GAS with Lucia), and the ASTA spin
!> density for SGUGA and Faroald.
subroutine Check_PDMs(iState,LSM,nVec,C,I12,n1,Rho1,n2,Rho2)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, CI_Chk_PDM_Faro, CI_Chk_PDM_SG, CI_Chk_SD_ASTA, CI_SG_Avail, &
                           CI_SG_Spin_Density_ASTA, CI_SG_Spin_PDM_G, CI_SG_Two_PDM_G
  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Input, only: iDoGAS
  use CI_Lib_Defs, only: CI_Lib_Selected, CI_Lib_Lucia, CI_Lib_SGUGA
  use Lucia_Interface, only: Lucia_Util
  use Lucia_Data, only: DSTmp, DTmp, PATmp, PTmp
  use sguga, only: CIS, EXS, SGS
  use SGUGA_Densities, only: sg_one_pdm
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, LSM, nVec, I12, n1, n2
  real(kind=wp), intent(in) :: C(nVec), Rho1(n1), Rho2(n2)
  integer(kind=iwp) :: NAC, NACPAR, NACPR2
  real(kind=wp) :: Dev(4)
  logical(kind=iwp) :: Done
  real(kind=wp), allocatable :: CIL(:), DS(:), PAS(:), PS(:), R2(:), SD(:), SDA(:), SDS(:)
  real(kind=wp), parameter :: ThrD = 1.0e-10_wp

  NAC = CI_Lib_Space(iState)%nLev
  NACPAR = nTri_Elem(NAC)
  NACPR2 = nTri_Elem(NACPAR)
  call mma_allocate(SD,max(1,NACPAR),Label='SD')
  call mma_allocate(R2,max(1,n2),Label='R2')
  call Raw_Dens(iState,1,1,LSM,LSM,nVec,C,nVec,C,size(SD),SD,size(R2),R2,Packed=.true.)
  call mma_deallocate(R2)
  call mma_allocate(CIL,nVec,Label='CIL')
  if (CI_Lib_Selected /= CI_Lib_Lucia) then
    call CI_Basis_To_Lucia(iState,LSM,nVec,C,CIL)
    call mma_allocate(DTMP,NAC**2,Label='DTMP')
    call mma_allocate(DSTMP,NAC**2,Label='DSTMP')
    call mma_allocate(PTMP,NACPR2,Label='PTMP')
    call mma_allocate(PATMP,NACPR2,Label='PATMP')
    call Lucia_Util('Densi',CI_Vector=CIL)
    ! (only the triangles are defined; no spin density for singlets)
    Dev(:) = Zero
    Dev(1) = maxval(abs(DTMP(1:NACPAR)-Rho1(1:NACPAR)))
    if (CI_Lib_Space(iState)%iSpin > 1) Dev(2) = maxval(abs(DSTMP(1:NACPAR)-SD(1:NACPAR)))
    if (I12 == 2) then
      Dev(3) = maxval(abs(PTMP(1:NACPR2)-Rho2(1:NACPR2)))
      Dev(4) = maxval(abs(PATMP(1:NACPR2)-Rho2(NACPR2+1:2*NACPR2)))
    end if
    Done = CI_Chk_PDM_SG
    call CI_Chk_Report_Dev('CI_Lib_Dens: D, SD, P, PA vs Lucia',maxval(Dev),ThrD,Done)
    call mma_deallocate(PATMP)
    call mma_deallocate(PTMP)
    call mma_deallocate(DSTMP)
    call mma_deallocate(DTMP)
  end if
  if ((CI_Lib_Selected /= CI_Lib_SGUGA) .and. CI_SG_Avail() .and. ((CI_Lib_Selected /= CI_Lib_Lucia) .or. (.not. iDoGAS))) then
    call CI_Basis_To_File(iState,LSM,nVec,C,CIL)
    call mma_allocate(DS,NACPAR,Label='DS')
    call mma_allocate(SDS,NACPAR,Label='SDS')
    call mma_allocate(PS,NACPR2,Label='PS')
    call mma_allocate(PAS,NACPR2,Label='PAS')
    call sg_one_pdm(SGS(iState),CIS(iState),EXS(iState),CIL,nVec,LSM,DS,NACPAR)
    Dev(:) = Zero
    Dev(1) = maxval(abs(DS(:)-Rho1(1:NACPAR)))
    if (CI_Lib_Selected /= CI_Lib_Lucia) then
      call CI_SG_Spin_PDM_G(iState,CIL,nVec,LSM,DS,SDS,NACPAR)
      Dev(2) = maxval(abs(SDS(:)-SD(1:NACPAR)))
    end if
    if (I12 == 2) then
      call CI_SG_Two_PDM_G(iState,CIL,nVec,LSM,PS,PAS,NACPR2)
      Dev(3) = maxval(abs(PS(1:NACPR2)-Rho2(1:NACPR2)))
      Dev(4) = maxval(abs(PAS(1:NACPR2)-Rho2(NACPR2+1:2*NACPR2)))
    end if
    Done = CI_Chk_PDM_Faro
    call CI_Chk_Report_Dev('CI_Lib_Dens: D, SD, P, PA vs SGUGA',maxval(Dev),ThrD,Done)
    call mma_deallocate(PAS)
    call mma_deallocate(PS)
    call mma_deallocate(SDS)
    call mma_deallocate(DS)
  end if
  if ((CI_Lib_Selected /= CI_Lib_Lucia) .and. CI_SG_Avail()) then
    ! the ASTA spin density (split-graph order)
    call CI_Basis_To_File(iState,LSM,nVec,C,CIL)
    call mma_allocate(SDA,NACPAR,Label='SDA')
    call CI_SG_Spin_Density_ASTA(iState,LSM,nVec,CIL,NACPAR,SDA)
    Done = CI_Chk_SD_ASTA
    call CI_Chk_Report_Dev('CI_Lib_Dens: spin density, ASTA T_pq vs the CI library',maxval(abs(SDA(:)-SD(1:NACPAR))),ThrD,Done)
    CI_Chk_SD_ASTA = .true.
    call mma_deallocate(SDA)
  end if
  call mma_deallocate(CIL)
  call mma_deallocate(SD)

end subroutine Check_PDMs

!> CI development (SIGMA): the one-body transition density Rho1 (spin-free,
!> layout of CI_Lib_Dens) of Bra and Ket and its spin part (computed here)
!> against Lucia (no spin density for singlets) and SGUGA (sg_trans_pdm,
!> whose layout is the transpose), those not in use; in the level layout
!> T(p+(q-1)*nAct) = <Bra|E_qp|Ket> (Lucia's).
subroutine Check_T1DM(iState,LSM,nVec,Bra_Vec,Ket_Vec,n1,Rho1,n2)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, CI_SG_Avail, CI_SG_Trans_PDM_G, CI_Chk_T1DM_Faro, CI_Chk_T1DM_SG
  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Defs, only: CI_Lib_Selected, CI_Lib_Lucia, CI_Lib_SGUGA
  use Lucia_Interface, only: Lucia_Util
  use Lucia_Data, only: DSTmp, DTmp
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, LSM, nVec, n1, n2
  real(kind=wp), intent(in) :: Bra_Vec(nVec), Ket_Vec(nVec), Rho1(n1)
  integer(kind=iwp) :: NAC, p, q
  real(kind=wp) :: Dev(2)
  logical(kind=iwp) :: Done
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: BL(:), KL(:), R1(:), R2(:), T1DM(:), T1SDM(:), TS(:), TSS(:)
  real(kind=wp), parameter :: ThrD = 1.0e-10_wp

  NAC = CI_Lib_Space(iState)%nLev
  call mma_allocate(L2O,max(1,NAC),Label='L2O')
  call CI_Lib_Level_Map(iState,NAC,L2O)
  call mma_allocate(T1DM,max(1,NAC**2),Label='T1DM')
  call mma_allocate(T1SDM,max(1,NAC**2),Label='T1SDM')
  call mma_allocate(R1,max(1,NAC**2),Label='R1')
  call mma_allocate(R2,max(1,n2),Label='R2')
  call Raw_Dens(iState,1,1,LSM,LSM,nVec,Bra_Vec,nVec,Ket_Vec,size(R1),R1,size(R2),R2)
  do q=1,NAC
    do p=1,NAC
      T1DM(p+(q-1)*NAC) = Rho1(L2O(q)+(L2O(p)-1)*NAC)
      T1SDM(p+(q-1)*NAC) = R1(L2O(q)+(L2O(p)-1)*NAC)
    end do
  end do
  call mma_deallocate(R2)
  call mma_deallocate(R1)
  call mma_deallocate(L2O)

  call mma_allocate(BL,nVec,Label='BL')
  call mma_allocate(KL,nVec,Label='KL')
  if (CI_Lib_Selected /= CI_Lib_Lucia) then
    call CI_Basis_To_Lucia(iState,LSM,nVec,Bra_Vec,BL)
    call CI_Basis_To_Lucia(iState,LSM,nVec,Ket_Vec,KL)
    call mma_allocate(DTMP,NAC*NAC,Label='DTMP')
    call mma_allocate(DSTMP,NAC*NAC,Label='DSTMP')
    call Lucia_Util('Densi',CI_Vector=KL(:),RVec=BL(:))
    Dev(1) = maxval(abs(DTMP(1:NAC*NAC)-T1DM(1:NAC*NAC)))
    Dev(2) = Zero
    if (CI_Lib_Space(iState)%iSpin > 1) Dev(2) = maxval(abs(DSTMP(1:NAC*NAC)-T1SDM(1:NAC*NAC)))
    Done = CI_Chk_T1DM_Faro
    call CI_Chk_Report_Dev('CI_Lib_Dens: T1DM, T1SDM vs Lucia',maxval(Dev),ThrD,Done)
    call mma_deallocate(DSTMP)
    call mma_deallocate(DTMP)
  end if
  if ((CI_Lib_Selected /= CI_Lib_SGUGA) .and. CI_SG_Avail()) then
    call CI_Basis_To_File(iState,LSM,nVec,Bra_Vec,BL)
    call CI_Basis_To_File(iState,LSM,nVec,Ket_Vec,KL)
    call mma_allocate(TS,NAC*NAC,Label='TS')
    call mma_allocate(TSS,NAC*NAC,Label='TSS')
    call CI_SG_Trans_PDM_G(iState,BL,KL,nVec,LSM,NAC,TS,TSS,.true.)
    Dev(:) = Zero
    do q=1,NAC
      do p=1,NAC
        Dev(1) = max(Dev(1),abs(T1DM(p+(q-1)*NAC)-TS(q+(p-1)*NAC)))
        if (CI_Lib_Selected /= CI_Lib_Lucia) Dev(2) = max(Dev(2),abs(T1SDM(p+(q-1)*NAC)-TSS(q+(p-1)*NAC)))
      end do
    end do
    Done = CI_Chk_T1DM_SG
    call CI_Chk_Report_Dev('CI_Lib_Dens: T1DM, T1SDM vs SGUGA',maxval(Dev),ThrD,Done)
    call mma_deallocate(TSS)
    call mma_deallocate(TS)
  end if
  CI_Chk_T1DM_Faro = .true.
  CI_Chk_T1DM_SG = .true.
  call mma_deallocate(KL)
  call mma_deallocate(BL)
  call mma_deallocate(T1SDM)
  call mma_deallocate(T1DM)

end subroutine Check_T1DM

!> The diagonal of H + ECore of the core (CI_Lib_Diag of CI_Lib_Core).
!> For the configuration average of the state from the CI input (the
!> diagonal of the CI solver, no core energy) the development options of
!> MOLCAS_CI_VERIFY: SIGMA compares it with the other libraries and with the
!> exact SGUGA diagonal (exactly for the single-determinant configurations;
!> with nPMax, the CSFs of the P space, also the selection of the P space by
!> either diagonal, reported, not enforced); HDIAG replaces it by the exact
!> SGUGA diagonal.
subroutine CI_Lib_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag,nPMax)

  use CI_Lib_Checks, only: CI_SG_Avail, CI_Tim_IHI, CI_Tim_Off, CI_Tim_On
  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Defs, only: CI_Diag_Cnf_Avg, CI_Lib_Selected, CI_Lib_Faroald, CI_Lib_Lucia, CI_Lib_SGUGA
  use CI_Verify_Options, only: CI_Verify
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, Mode, iSym, nTU, nTUVX, n
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp), intent(in), optional :: nPMax
  real(kind=wp) :: TimT(2)
  real(kind=wp), allocatable :: DL(:)

  call CI_Tim_On(TimT)
  call Raw_Diag(iState,Mode,iSym,nTU,TU,nTUVX,TUVX,ECore,n,Diag)
  if (CI_Lib_Input_State(iState) .and. (Mode == CI_Diag_Cnf_Avg) .and. (iSym == CI_Lib_Space(iState)%STSym) .and. &
      (ECore == Zero)) then
    select case (CI_Lib_Selected)
      case (CI_Lib_SGUGA)
        call Check_Diag_SG(iState,n,Diag,nTU,TU,nTUVX,TUVX)
      case (CI_Lib_Faroald)
        if (CI_Verify('SIGMA')) call Check_Diag_Faro(iState,n,Diag,nTU,TU,nTUVX,TUVX)
      case (CI_Lib_Lucia)
        if (CI_Verify('SIGMA') .and. CI_SG_Avail()) call Check_Diag_Lucia(iState,n,Diag,nTU,TU,nTUVX,TUVX,nPMax)
        if (CI_Verify('HDIAG') .and. CI_SG_Avail()) then
          ! the exact diagonal instead (in the Lucia order)
          call mma_allocate(DL,n,Label='DL')
          call Exact_Diag_Lucia(iState,n,nTU,TU,nTUVX,TUVX,DL)
          Diag(:) = DL(:)
          call mma_deallocate(DL)
          write(u6,'(1X,A)') 'CI_Lib_Diag: exact SGUGA diagonal used (MOLCAS_CI_VERIFY=HDIAG)'
        end if
    end select
  end if
  call CI_Tim_Off(CI_Tim_IHI,TimT)

end subroutine CI_Lib_Diag

!> CI development (SIGMA): Lucia's configuration-average diagonal of the CI
!> symmetry of state iState in the Lucia basis (CSDIAG of its determinant
!> diagonal)
subroutine Lucia_Avg_Diag(iState,n,nTU,TU,nTUVX,TUVX,DL)

  use CI_Lib_Info, only: CI_Lib_Space
  use Lucia_Interface, only: Lucia_Util
  use lucia_data, only: NCNFTP, NCSFTP, NDTFTP, NTYP, SDREO
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: DL(n)
  integer(kind=iwp) :: nDet
  real(kind=wp), allocatable :: DD(:)
  integer(kind=iwp), parameter :: IPRINT = 0

  ! the determinants of Lucia
  nDet = maxval(CI_Lib_Space(iState)%nDet(:))
  call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
  call mma_allocate(DD,nDet,label='DETDIA')
  call get_diag(DD,nDet)
  call CSDIAG(n,nDet,DL,DD,NCNFTP(1,CI_Lib_Space(iState)%STSym),NTYP,SDREO,NDTFTP,NCSFTP,IPRINT)
  call mma_deallocate(DD)

end subroutine Lucia_Avg_Diag

!> A diagonal in the SGUGA CSF order to the Lucia order: a permutation,
!> without the CSF phases (CI_SG_ReOrd of the index vector: Idx(ICSF) = +-ISG).
subroutine Diag_To_Lucia(iState,n,DSG,DL)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_SG_Lucia_Bridge, only: CI_SG_ReOrd
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, n
  real(kind=wp), intent(in) :: DSG(n)
  real(kind=wp), intent(out) :: DL(n)
  integer(kind=iwp) :: j
  real(kind=wp), allocatable :: Idx(:), Map(:)

  call mma_allocate(Map,n,Label='Map')
  call mma_allocate(Idx,n,Label='Idx')
  do j=1,n
    Map(j) = real(j,kind=wp)
  end do
  call CI_SG_ReOrd(iState,CI_Lib_Space(iState)%STSym,1,n,Map,Idx)
  do j=1,n
    DL(j) = DSG(nint(abs(Idx(j))))
  end do
  call mma_deallocate(Idx)
  call mma_deallocate(Map)

end subroutine Diag_To_Lucia

!> The exact diagonal <I|H|I> (sg_hdiag_fast) in the Lucia CSF order.
subroutine Exact_Diag_Lucia(iState,n,nTU,TU,nTUVX,TUVX,DL)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Checks, only: CI_SG_HDiag_G
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: DL(n)
  real(kind=wp), allocatable :: DSG(:)

  call mma_allocate(DSG,n,Label='DSG')
  call CI_SG_HDiag_G('FAST',iState,n,CI_Lib_Space(iState)%STSym,TUVX,nTUVX,TU,nTU,DSG)
  call CI_HEXS_Proj(CI_HEXS_Keep_SG,DSG)
  call Diag_To_Lucia(iState,n,DSG,DL)
  call mma_deallocate(DSG)

end subroutine Exact_Diag_Lucia

!> SGCI: the configuration average of the core (SG_Cnf_Avg_Diag, so that
!> SGCI follows the iterations of Lucia), or with MOLCAS_CI_VERIFY=HDIAG the
!> exact <I|H|I> (sg_hdiag_fast). Under SIGMA the average is checked against
!> Lucia and the exact diagonal against sg_hdiag_ref (small spaces).
subroutine Check_Diag_SG(iState,n,Diag,nTU,TU,nTUVX,TUVX)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Checks, only: CI_Chk_IHI_Avg, CI_Chk_IHI_SG, CI_SG_HDiag_G
  use CI_Lib_Defs, only: CI_Diag_Exact
  use CI_Verify_Options, only: CI_Verify
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(inout) :: Diag(n)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  integer(kind=iwp), parameter :: MxRef = 500
  real(kind=wp), parameter :: ThrDia = 1.0e-10_wp
  real(kind=wp) :: MaxDev
  real(kind=wp), allocatable :: DAvgL(:), DEx(:), DL(:), DRef(:)

  if (CI_Verify('SIGMA')) then
    ! Lucia's configuration-average diagonal, in the Lucia basis
    call mma_allocate(DL,n,Label='DL')
    call Lucia_Avg_Diag(iState,n,nTU,TU,nTUVX,TUVX,DL)
    call mma_allocate(DAvgL,n,Label='DAvgL')
    call Diag_To_Lucia(iState,n,Diag,DAvgL)
    MaxDev = maxval(abs(DAvgL(:)-DL(:)))
    if (.not. CI_Chk_IHI_Avg) &
      write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Diag (SGCI): SGUGA configuration-average diagonal verified against Lucia, ',n, &
                                     ' CSFs, max. error ',MaxDev
    CI_Chk_IHI_Avg = .true.
    if (MaxDev > ThrDia*max(One,maxval(abs(DL(:))))) then
      write(u6,*) 'CI_Lib_Diag: the SGUGA configuration-average diagonal differs from Lucia'
      call Abend()
    end if
    call mma_deallocate(DAvgL)
    call mma_deallocate(DL)
  end if
  if (.not. (CI_Verify('HDIAG') .or. CI_Verify('SIGMA'))) return
  ! the exact diagonal (HDIAG, and its check under SIGMA)
  call mma_allocate(DEx,n,Label='DExact')
  call Raw_Diag(iState,CI_Diag_Exact,CI_Lib_Space(iState)%STSym,nTU,TU,nTUVX,TUVX,Zero,n,DEx)
  if (CI_Verify('HDIAG')) then
    Diag(:) = DEx(:)
    write(u6,'(1X,A)') 'CI_Lib_Diag (SGCI): exact SGUGA diagonal used (MOLCAS_CI_VERIFY=HDIAG)'
  end if
  if (CI_Verify('SIGMA') .and. (n <= MxRef)) then
    call mma_allocate(DRef,n,Label='DRef')
    call CI_SG_HDiag_G('REF',iState,n,CI_Lib_Space(iState)%STSym,TUVX,nTUVX,TU,nTU,DRef)
    call CI_HEXS_Proj(CI_HEXS_Keep_SG,DRef)
    MaxDev = maxval(abs(DEx(:)-DRef(:)))
    if (.not. CI_Chk_IHI_SG) &
      write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Diag (SGCI): sg_hdiag_fast verified against sg_hdiag_ref, ',n,' CSFs, max. error ', &
                                     MaxDev
    CI_Chk_IHI_SG = .true.
    if (MaxDev > ThrDia*max(One,maxval(abs(DRef(:))))) then
      write(u6,*) 'CI_Lib_Diag: sg_hdiag_fast differs from sg_hdiag_ref'
      call Abend()
    end if
    call mma_deallocate(DRef)
  end if
  call mma_deallocate(DEx)

end subroutine Check_Diag_SG

!> FARO (SIGMA): the configuration average of Faroald (from the
!> occupations) checked once against the average of CI_Faro_HDet over the
!> determinants (M_S = S) of each kept configuration, and against Lucia and
!> SGUGA.
subroutine Check_Diag_Faro(iState,n,Diag,nTU,TU,nTUVX,TUVX)

  use CI_Lib_Checks, only: CI_Basis_To_Lucia, CI_Chk_IHI_Avg, CI_SG_Avail
  use CI_Lib_Info, only: CI_Lib_Space
  use Faroald_CITrans, only: Faro_conf_dets
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(in) :: Diag(n), TU(nTU), TUVX(nTUVX)
  real(kind=wp), parameter :: ThrDia = 1.0e-10_wp
  integer(kind=iwp) :: i, k, k1, nc, nd
  real(kind=wp) :: E, MaxDev
  logical(kind=iwp), allocatable :: Kept(:)
  integer(kind=iwp), allocatable :: da(:), db(:)
  real(kind=wp), allocatable :: cf(:,:), DL(:), Idx(:), Map(:)

  if (.not. CI_Chk_IHI_Avg) then
    ! the average of the determinant energies, as reference (the removed
    ! configurations of HEXS/DEXS, zero in Diag, skipped)
    call CI_Faro_Ints_Act(1,nTU,TU,nTUVX,TUVX,.true.)
    call mma_allocate(Kept,n,Label='Kept')
    Kept(:) = .true.
    if (allocated(CI_HEXS_Keep)) then
      if (size(CI_HEXS_Keep) == n) Kept(:) = CI_HEXS_Keep(:)
    end if
    MaxDev = Zero
    k = 1
    do while (k <= n)
      call Faro_conf_dets(k,k1,nc,nd,da,db,cf)
      E = Zero
      do i=1,nd
        E = E+CI_Faro_HDet(da(i),db(i),da(i),db(i))
      end do
      if (Kept(k1)) MaxDev = max(MaxDev,abs(Diag(k1)-E/real(max(nd,1),kind=wp)))
      k = k1+nc
      call mma_deallocate(da)
      call mma_deallocate(db)
      call mma_deallocate(cf)
    end do
    call mma_deallocate(Kept)
    write(u6,'(1X,A,ES10.3)') 'CI_Lib_Diag (Faroald): occupation formula vs determinant average, max. difference ',MaxDev
    if (MaxDev > ThrDia*max(One,maxval(abs(Diag(1:n))))) then
      write(u6,*) 'CI_Lib_Diag: the Faroald occupation formula differs from the determinant average'
      call Abend()
    end if
  end if
  call mma_allocate(DL,n,Label='DL')
  call Lucia_Avg_Diag(iState,n,nTU,TU,nTUVX,TUVX,DL)
  ! the permutation Faroald -> Lucia basis
  call mma_allocate(Idx,n,Label='Idx')
  call mma_allocate(Map,n,Label='Map')
  do i=1,n
    Idx(i) = real(i,kind=wp)
  end do
  call CI_Basis_To_Lucia(iState,CI_Lib_Space(iState)%STSym,n,Idx,Map)
  MaxDev = Zero
  do i=1,n
    MaxDev = max(MaxDev,abs(Diag(nint(abs(Map(i))))-DL(i)))
  end do
  if (.not. CI_Chk_IHI_Avg) &
    write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Diag (Faroald): configuration-average diagonal verified against Lucia, ',n, &
                                   ' CSFs, max. error ',MaxDev
  CI_Chk_IHI_Avg = .true.
  if (MaxDev > ThrDia*max(One,maxval(abs(DL(:))))) then
    write(u6,*) 'CI_Lib_Diag: the Faroald configuration-average diagonal differs from Lucia'
    call Abend()
  end if
  call mma_deallocate(Map)
  call mma_deallocate(Idx)
  call mma_deallocate(DL)
  if (CI_SG_Avail()) call Check_Diag_Faro_SG(iState,n,Diag,nTU,TU,nTUVX,TUVX)

end subroutine Check_Diag_Faro

!> FARO (SIGMA): the Faroald configuration-average diagonal against that of
!> SGUGA (SG_Cnf_Avg_Diag, split-graph order of the files).
subroutine Check_Diag_Faro_SG(iState,n,Diag,nTU,TU,nTUVX,TUVX)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Chk_IHI_Faro_SG
  use CI_Lib_Info, only: CI_Lib_Space
  use sguga, only: SG_Cnf_Avg_Diag, SG_Free_Cnf_List, SG_Mk_Cnf_List, SGCnfList
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(in) :: Diag(n), TU(nTU), TUVX(nTUVX)
  real(kind=wp), parameter :: ThrDia = 1.0e-10_wp
  type(SGCnfList) :: Cnf
  integer(kind=iwp) :: i
  real(kind=wp) :: MaxDev
  real(kind=wp), allocatable :: DSG(:), Idx(:), Map(:)

  call SG_Mk_Cnf_List(iState,CI_Lib_Space(iState)%STSym,Cnf)
  call mma_allocate(DSG,n,Label='DSG')
  call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,DSG)
  call CI_HEXS_Proj(CI_HEXS_Keep_SG,DSG)
  call mma_allocate(Idx,n,Label='Idx')
  call mma_allocate(Map,n,Label='Map')
  do i=1,n
    Idx(i) = real(i,kind=wp)
  end do
  call CI_Basis_To_File(iState,CI_Lib_Space(iState)%STSym,n,Idx,Map)
  MaxDev = Zero
  do i=1,n
    MaxDev = max(MaxDev,abs(Diag(nint(abs(Map(i))))-DSG(i)))
  end do
  call CI_Chk_Report_Dev('CI_Lib_Diag (FARO): Faroald configuration-average diagonal vs SGUGA',MaxDev, &
                  ThrDia*max(One,maxval(abs(DSG(:)))),CI_Chk_IHI_Faro_SG)
  call mma_deallocate(Map)
  call mma_deallocate(Idx)
  call mma_deallocate(DSG)
  call SG_Free_Cnf_List(Cnf)

end subroutine Check_Diag_Faro_SG

!> RASSCF's Lucia (SIGMA): its configuration-average diagonal Diag (Lucia
!> order) against the exact SGUGA diagonal (sg_hdiag_fast, for small spaces
!> checked against sg_hdiag_ref): exactly for the single-determinant
!> configurations, where the average is exact, and, with nPMax, by the
!> configurations that the selection of the P space would pick with either
!> diagonal (reported, not enforced: the average is an approximation).
subroutine Check_Diag_Lucia(iState,n,Diag,nTU,TU,nTUVX,TUVX,nPMax)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Checks, only: CI_SG_HDiag_G
  use lucia_data, only: MINOP, NCNFTP, NCSFTP, NTYP
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, n, nTU, nTUVX
  real(kind=wp), intent(in) :: Diag(n), TU(nTU), TUVX(nTUVX)
  integer(kind=iwp), intent(in), optional :: nPMax
  integer(kind=iwp), parameter :: MxRef = 500, MxPrt = 10
  real(kind=wp), parameter :: ThrDia = 1.0e-10_wp
  integer(kind=iwp) :: i, iC, iCnf, iRk(3), iRkMax, iRkMin, iTyp, k, LSM, nBoth, nC, nCnf, nPrt, nSD, nSDErr, nSelL, nSelMean, &
                       nSelMin
  real(kind=wp) :: Dev, MaxDev
  integer(kind=iwp), allocatable :: nCSFCnf(:), nOpCnf(:)
  real(kind=wp), allocatable :: DEx(:), DRef(:), DSG(:), KeyL(:), KeyMean(:), KeyMin(:)
  logical(kind=iwp), allocatable :: SelL(:), SelMean(:), SelMin(:)

  LSM = CI_Lib_Space(iState)%STSym
  ! exact diagonal in the SGUGA order
  call mma_allocate(DSG,n,Label='DSG')
  call CI_SG_HDiag_G('FAST',iState,n,LSM,TUVX,nTUVX,TU,nTU,DSG)
  if (n <= MxRef) then
    call mma_allocate(DRef,n,Label='DRef')
    call CI_SG_HDiag_G('REF',iState,n,LSM,TUVX,nTUVX,TU,nTU,DRef)
    MaxDev = Zero
    do i=1,n
      Dev = abs(DSG(i)-DRef(i))
      MaxDev = max(MaxDev,Dev)
      if (Dev > ThrDia*max(One,abs(DRef(i)))) then
        write(u6,*) 'CI_Lib_Diag: sg_hdiag_fast differs from sg_hdiag_ref for CSF',i,DSG(i),DRef(i)
        call Abend()
      end if
    end do
    write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Diag: sg_hdiag_fast verified against sg_hdiag_ref, ',n,' CSFs, max. error ',MaxDev
    call mma_deallocate(DRef)
  end if

  ! to the Lucia CSF order (with HEXS/DEXS zero for the removed CSFs, as
  ! Lucia's)
  call CI_HEXS_Proj(CI_HEXS_Keep_SG,DSG)
  call mma_allocate(DEx,n,Label='DEx')
  call Diag_To_Lucia(iState,n,DSG,DEx)
  call mma_deallocate(DSG)

  ! configuration layout with the number of open shells (NDTFTP counts
  ! spin combinations, not determinants, for MS = 0)
  nCnf = sum(NCNFTP(1:NTYP,LSM))
  call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
  call mma_allocate(nOpCnf,nCnf,Label='nOpCnf')
  k = 0
  do iTyp=1,NTYP
    do iCnf=1,NCNFTP(iTyp,LSM)
      k = k+1
      nCSFCnf(k) = NCSFTP(iTyp)
      nOpCnf(k) = iTyp-1+MINOP
    end do
  end do

  ! keys per configuration: Lucia (the average, as the P-space selection
  ! uses it) and the exact lowest and mean CSF energies; exact check for
  ! single determinants: closed shells, or all open shells parallel (2S open
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
    KeyL(k) = Diag(iC+1)
    KeyMin(k) = minval(DEx(iC+1:iC+nC))
    KeyMean(k) = sum(DEx(iC+1:iC+nC))/real(nC,kind=wp)
    if (nOpCnf(k) == CI_Lib_Space(iState)%iSpin-1) then
      nSD = nSD+1
      Dev = abs(Diag(iC+1)-DEx(iC+1))
      MaxDev = max(MaxDev,Dev)
      if (Dev > ThrDia*max(One,abs(DEx(iC+1)))) then
        nSDErr = nSDErr+1
        if (nSDErr <= MxPrt) write(u6,*) 'CI_Lib_Diag: single-det. configuration',k,' CI solver:',Diag(iC+1),' exact:',DEx(iC+1)
      end if
    end if
    iC = iC+nC
  end do
  if (nSDErr > 0) then
    write(u6,*) 'CI_Lib_Diag: diagonal of single-determinant configurations differs from SGUGA'
    call Abend()
  end if
  write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Diag: diagonal verified against SGUGA for ',nSD,' single-determinant configurations, &
                                 &max. error ',MaxDev

  if (present(nPMax)) then
    ! the selection of the P space with each key
    call mma_allocate(SelL,nCnf,Label='SelL')
    call mma_allocate(SelMin,nCnf,Label='SelMin')
    call mma_allocate(SelMean,nCnf,Label='SelMean')
    call Sel_PHP(nCnf,KeyL,nCSFCnf,nPMax,SelL)
    call Sel_PHP(nCnf,KeyMin,nCSFCnf,nPMax,SelMin)
    call Sel_PHP(nCnf,KeyMean,nCSFCnf,nPMax,SelMean)
    nSelL = count(SelL)
    nSelMin = count(SelMin)
    nSelMean = count(SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A)') 'CI_Lib_Diag: explicit subspace (NSEL = ',nPMax,' of ',n,' CSFs, ',nCnf, &
                                           ' configurations) from the CI-solver diagonal: ',nSelL,' configurations'
    if (nPMax >= n) write(u6,'(1X,A)') '        (the whole space, the comparison is trivial)'
    nBoth = count(SelL .and. SelMin)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact lowest CSF energy per configuration: ',nSelMin,' configurations, ',nBoth, &
                                      ' in common, ',count(SelL .neqv. SelMin),' differ'
    nBoth = count(SelL .and. SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact mean CSF energy per configuration:   ',nSelMean,' configurations, ',nBoth, &
                                      ' in common, ',count(SelL .neqv. SelMean),' differ'
    ! The selection stops at the first configuration that does not fit into
    ! nPMax CSFs, so near the cutoff the sizes of the configurations decide;
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
        iRk(1) = CI_Cnf_Rank(nCnf,KeyL,k)
        iRk(2) = CI_Cnf_Rank(nCnf,KeyMin,k)
        iRk(3) = CI_Cnf_Rank(nCnf,KeyMean,k)
        iRkMin = min(iRkMin,minval(iRk))
        iRkMax = max(iRkMax,maxval(iRk))
        nPrt = nPrt+1
        if (nPrt <= MxPrt) write(u6,'(9X,I8,3L2,2I6,3I8)') k,SelL(k),SelMin(k),SelMean(k),nCSFCnf(k),nOpCnf(k),iRk(:)
      end do
      write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        differing configurations are ranked ',iRkMin,' to ',iRkMax,' (',nSelL, &
                                        ' selected by the CI solver)'
    end if
    call mma_deallocate(SelMean)
    call mma_deallocate(SelMin)
    call mma_deallocate(SelL)
  end if

  call mma_deallocate(KeyMean)
  call mma_deallocate(KeyMin)
  call mma_deallocate(KeyL)
  call mma_deallocate(nOpCnf)
  call mma_deallocate(nCSFCnf)
  call mma_deallocate(DEx)

end subroutine Check_Diag_Lucia

!> Position of Key(k) in the order of increasing Key (ties by index).
function CI_Cnf_Rank(nKey,Key,k) result(iRank)

  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: nKey, k
  real(kind=wp), intent(in) :: Key(nKey)
  integer(kind=iwp) :: iRank, l

  iRank = 1
  do l=1,nKey
    if ((Key(l) < Key(k)) .or. ((Key(l) == Key(k)) .and. (l < k))) iRank = iRank+1
  end do

end function CI_Cnf_Rank

!> The configurations (nCSFCnf CSFs each) of a P space of at most nPMax
!> CSFs with the configuration energies Key, as RASSCF's Lucia selects
!> them: the lowest first (ties within Acc in the order of the
!> configurations), until the next one does not fit.
subroutine Sel_PHP(nKey,Key,nCSFCnf,nPMax,Sel)

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: nKey, nCSFCnf(nKey), nPMax
  real(kind=wp), intent(in) :: Key(nKey)
  logical(kind=iwp), intent(out) :: Sel(nKey)
  integer(kind=iwp) :: k, kMin, nCMin, nP, nPCnf
  real(kind=wp) :: XMax, XMin
  real(kind=wp), allocatable :: Wrk(:)
  real(kind=wp), parameter :: Acc = 1.0e-13_wp

  call mma_allocate(Wrk,nKey,Label='Wrk')
  Wrk(:) = Key(:)
  XMax = maxval(abs(Key))
  Sel(:) = .false.
  nP = 0
  nPCnf = 0
  do
    XMin = XMax+One
    kMin = 0
    nCMin = 0
    do k=1,nKey
      if (Wrk(k)+Acc < XMin) then
        XMin = Wrk(k)
        kMin = k
        nCMin = nCSFCnf(k)
      end if
    end do
    if ((kMin == 0) .or. (nP+nCMin > nPMax)) exit
    Sel(kMin) = .true.
    nP = nP+nCMin
    nPCnf = nPCnf+1
    Wrk(kMin) = XMax+One
    if (nPCnf >= nKey) exit
  end do
  call mma_deallocate(Wrk)

end subroutine Sel_PHP

!> S = H C of the core (CI_Lib_Sigma of CI_Lib_Core, see there). Under
!> MOLCAS_CI_VERIFY=SIGMA, for H in the symmetry of the state from the CI
!> input, compared with Lucia (SGCI, FARO) and with SGUGA, whose sigma
!> vector is then used (with RAS partitions from sg_h_psi_ras; for CAS from
!> sg_h_psi, checked once against sg_h_psi_ras).
subroutine CI_Lib_Sigma(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Lib_Checks, only: CI_Tim_HPsi, CI_Tim_Off, CI_Tim_On
  use CI_Verify_Options, only: CI_Verify
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, iOpSpin, iCSym, iSSym, nTU, nTUVX, nC, nS
  logical(kind=iwp), intent(in) :: Square, Have_2_el
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), C(nC)
  real(kind=wp), intent(out) :: S(nS)
  real(kind=wp) :: TimT(2)

  call CI_Tim_On(TimT)
  call Raw_Sigma(iState,iOpSpin,iCSym,iSSym,Square,Have_2_el,nTU,TU,nTUVX,TUVX,nC,C,nS,S)
  if (CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (iOpSpin == 0) .and. (.not. Square) .and. Have_2_el .and. &
      (iCSym == CI_Lib_Space(iState)%STSym) .and. (iSSym == iCSym) .and. (nS == nC)) &
    call Check_Sigma(iState,iCSym,nC,C,S,nTU,TU,nTUVX,TUVX)
  call CI_Tim_Off(CI_Tim_HPsi,TimT)

end subroutine CI_Lib_Sigma

!> The sigma vector in the determinant basis of RASSCF's Lucia (CASVB): the
!> core's CI_Lib_Sigma_Det (see there). Under MOLCAS_CI_VERIFY=SIGMA
!> compared with the SGUGA sigma vector, which is then used (CAS without
!> symmetry only: the SD/CSF codec exists for CAS, the comparison in the
!> determinant basis needs one symmetry).
subroutine CI_Lib_Sigma_Det(iState,iSym,nTU,TU,nTUVX,TUVX,nDet,C,S)

  use CI_Lib_Checks, only: CI_SG_HPsi, CI_Tim_HPsi, CI_Tim_Off, CI_Tim_On
  use CI_Verify_Options, only: CI_Verify
  use CI_Lib_Info, only: CI_Lib_Space
  use CI_SD_CSF_Codec, only: CSF_to_SD, SD_CSF_Codec, SD_to_CSF
  use CI_SG_Lucia_Bridge, only: CI_SG_ReOrd
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nDet
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout) :: C(nDet)
  real(kind=wp), intent(out) :: S(nDet)
  integer(kind=iwp) :: iDet, nCSF
  real(kind=wp) :: TimT(2)
  real(kind=wp), allocatable :: T1(:), T2(:)

  call CI_Tim_On(TimT)
  call Raw_Sigma_Det(iState,iSym,nTU,TU,nTUVX,TUVX,nDet,C,S)
  if (CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (CI_Lib_Space(iState)%nRsPrt == 1) .and. &
      (CI_Lib_Space(iState)%nSym == 1)) then
    nCSF = CI_Lib_Space(iState)%nCSF(iSym)
    call mma_allocate(T1,nDet,Label='T1')
    call mma_allocate(T2,nDet,Label='T2')
    ! Lucia determinants -> Lucia CSFs -> split-graph CSFs
    T2(:) = C(:)
    T1(:) = Zero
    call SD_CSF_Codec(T2,T1,SD_to_CSF,iSym)
    call CI_SG_ReOrd(iState,iSym,0,nCSF,T1,T2)
    call CI_SG_HPsi(iState,T2,nCSF,iSym,T1,TUVX,nTUVX,TU,nTU,.false.)
    ! and back
    call CI_SG_ReOrd(iState,iSym,1,nCSF,T1,T2)
    call SD_CSF_Codec(T1,T2,CSF_to_SD,iSym)
    do iDet=1,nDet
      if (abs(T1(iDet)-S(iDet)) > 1.0e-10_wp) then
        write(u6,*) 'CI_Lib_Sigma_Det: Abs(temp1(iDet)-Sigma_vec(iDet))>1.0e-10_wp)'
        write(u6,*) ' iDet=',iDet
        write(u6,*) 'temp1(iDet)=    ',T1(iDet)
        write(u6,*) 'Sigma_vec(iDet)=',S(iDet)
        call Abend()
      end if
    end do
    S(:) = T1(:)
    call mma_deallocate(T2)
    call mma_deallocate(T1)
  end if
  call CI_Tim_Off(CI_Tim_HPsi,TimT)

end subroutine CI_Lib_Sigma_Det

!> CI development (SIGMA): the sigma vector S of C (H of the state iState
!> in its symmetry LSM) against Lucia (SGCI, FARO) and SGUGA; S is replaced
!> by SGUGA's (not for SGCI), projected for HEXS/DEXS.
subroutine Check_Sigma(iState,LSM,n,C,S,nTU,TU,nTUVX,TUVX)

  use CI_Lib_Info, only: CI_Lib_Space
  use CI_Input, only: iDoGAS
  use CI_Lib_Checks, only: CI_Basis_To_File, CI_File_To_Basis, CI_Faro_CI, CI_Chk_HPsi_Faro, CI_Chk_HPsi_RAS, &
                           CI_Chk_HPsi_SG, CI_Chk_Ref_Sigma, CI_SG_Avail, CI_SG_CI, CI_SG_HPsi
  use CI_SG_Lucia_Bridge, only: CI_SG_ReOrd
  use sguga, only: CIS, EXS, SGS
  use SGUGA_Sigma, only: sg_h_psi, sg_h_psi_ras
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, LSM, n, nTU, nTUVX
  real(kind=wp), intent(in) :: C(n), TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout) :: S(n)
  real(kind=wp) :: Check_H, Check_HRef, DMax
  logical(kind=iwp) :: Lucia_Done, SG_Done
  real(kind=wp), allocatable :: CM(:), SG_PSI(:), SG_RAS(:), SG_SGM(:)

  ! the vector the library got (HEXS/DEXS: projected)
  call mma_allocate(CM,n,Label='CM')
  CM(:) = C(:)
  call CI_HEXS_Proj(CI_HEXS_Keep,CM)

  ! SGCI: the SGUGA sigma vector against that of Lucia (vector reordered to
  ! the Lucia basis and the result back)
  if (CI_SG_CI()) then
    call mma_allocate(SG_PSI,max(n,CI_Lib_Space(iState)%nDet(LSM)),Label='SG_PSI')
    call mma_allocate(SG_SGM,max(n,CI_Lib_Space(iState)%nDet(LSM)),Label='SG_SGM')
    call CI_SG_ReOrd(iState,LSM,1,n,CM,SG_PSI)
    call Lucia_Sigma_Ref(iState,LSM,n,SG_PSI,SG_SGM,nTU,TU,nTUVX,TUVX)
    call CI_SG_ReOrd(iState,LSM,0,n,SG_SGM,SG_PSI)
    DMax = maxval(abs(SG_PSI(1:n)-S(:)))
    if (.not. CI_Chk_HPsi_SG) then
      write(u6,'(1X,A,ES10.3)') 'CI_Lib_Sigma (SGCI): SGUGA sigma vs Lucia, max. abs. difference',DMax
      CI_Chk_HPsi_SG = .true.
    end if
    if (DMax > 1.0e-10_wp*max(One,maxval(abs(S(:))))) then
      write(u6,'(1X,A,ES10.3)') 'SGUGA error in H|Psi> (SGCI), max. abs. difference',DMax
      call Abend()
    end if
    call mma_deallocate(SG_SGM)
    call mma_deallocate(SG_PSI)
    call mma_deallocate(CM)
    return
  end if

  ! FARO: the Faroald sigma vector against Lucia (and below against SGUGA)
  if (CI_Faro_CI()) then
    call mma_allocate(SG_SGM,n,Label='Lucia_SGM')
    call CI_Chk_Ref_Sigma('Lucia',iState,LSM,n,CM,SG_SGM,nTU,TU,nTUVX,TUVX)
    DMax = maxval(abs(SG_SGM(:)-S(:)))
    Lucia_Done = CI_Chk_HPsi_Faro
    call CI_Chk_Report_Dev('CI_Lib_Sigma (FARO): Faroald sigma vs Lucia',DMax,1.0e-10_wp*max(One,maxval(abs(S(:)))),Lucia_Done)
    call mma_deallocate(SG_SGM)
  end if

  ! the SGUGA sigma vector, which is then used: with RAS partitions from
  ! sg_h_psi_ras (no intermediate states outside the RAS space); for CAS from
  ! sg_h_psi, checked once per run against sg_h_psi_ras. With HEXS/DEXS
  ! projected (P H P, CI_HEXS_Keep_SG).
  if (CI_SG_Avail()) then
    call mma_allocate(SG_PSI,n,Label='SG_PSI')
    call CI_Basis_To_File(iState,LSM,n,CM,SG_PSI)
    call CI_HEXS_Proj(CI_HEXS_Keep_SG,SG_PSI)
    Check_HRef = CheckSum(S,n)
    call mma_allocate(SG_SGM,n,Label='SG_SGM')
    if ((CI_Lib_Space(iState)%nRsPrt == 1) .and. (.not. iDoGAS)) then
      call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),SG_PSI,n,LSM,SG_SGM,TUVX,nTUVX,TU,nTU)
      if (.not. CI_Chk_HPsi_RAS) then
        call mma_allocate(SG_RAS,n,Label='SG_RAS')
        call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),SG_PSI,n,LSM,SG_RAS,TUVX,nTUVX,TU,nTU)
        DMax = maxval(abs(SG_RAS(:)-SG_SGM(:)))
        write(u6,'(1X,A,ES10.3)') 'CI_Lib_Sigma: sg_h_psi_ras vs sg_h_psi, max. abs. difference',DMax
        if (DMax > 1.0e-12_wp*max(One,maxval(abs(SG_SGM(:))))) then
          write(u6,*) 'SGUGA error in sg_h_psi_ras'
          call Abend()
        end if
        call mma_deallocate(SG_RAS)
        CI_Chk_HPsi_RAS = .true.
      end if
    else
      call CI_SG_HPsi(iState,SG_PSI,n,LSM,SG_SGM,TUVX,nTUVX,TU,nTU,.true.)
    end if
    call CI_HEXS_Proj(CI_HEXS_Keep_SG,SG_SGM)
    call CI_File_To_Basis(iState,LSM,n,SG_SGM,SG_PSI)
    if ((CI_Lib_Space(iState)%nRsPrt /= 1) .or. iDoGAS) then
      DMax = maxval(abs(SG_PSI(:)-S(:)))
      if (.not. CI_Chk_HPsi_RAS) &
        write(u6,'(1X,A,ES10.3)') 'CI_Lib_Sigma: sg_h_psi_ras vs CI back end (RAS), max. abs. difference',DMax
      CI_Chk_HPsi_RAS = .true.
      if (DMax > 1.0e-10_wp*max(One,maxval(abs(S(:))))) then
        write(u6,'(1X,A,ES10.3)') 'SGUGA error in sg_h_psi_ras (RAS), max. abs. difference',DMax
        call Abend()
      end if
    end if
    if (CI_Faro_CI()) then
      SG_Done = CI_Chk_HPsi_Faro
      call CI_Chk_Report_Dev('CI_Lib_Sigma (FARO): Faroald sigma vs SGUGA',maxval(abs(SG_PSI(:)-S(:))), &
                      1.0e-10_wp*max(One,maxval(abs(S(:)))),SG_Done)
      CI_Chk_HPsi_Faro = .true.
    end if
    Check_H = CheckSum(SG_PSI,n)
    if (abs(Check_HRef-Check_H)/n > 1.0e-12_wp) then
      write(u6,*) 'SGUGA error in H|Psi>'
      write(u6,*) 'Check_HRef=',Check_HRef
      write(u6,*) 'Check_H   =',Check_H
      call RecPrt('Sigma_Vec',' ',S,1,n)
      call RecPrt('SG_SGM',' ',SG_PSI,1,n)
      call Abend()
    end if
    S(1:n) = SG_PSI(1:n)
    call CI_HEXS_Proj(CI_HEXS_Keep,S)
    call mma_deallocate(SG_SGM)
    call mma_deallocate(SG_PSI)
  end if
  call mma_deallocate(CM)

end subroutine Check_Sigma

!> Lucia's sigma vector Sout of the CSF vector Cin (Lucia basis, symmetry
!> LSM of state iState), the reference of SGCI; Cin and Sout hold
!> max(n,nDet) elements (CSDTVC).
subroutine Lucia_Sigma_Ref(iState,LSM,n,Cin,Sout,nTU,TU,nTUVX,TUVX)

  use CI_Lib_Info, only: CI_Lib_Space
  use Lucia_Interface, only: Lucia_Util
  use lucia_data, only: IREFSM, Sigma_on_disk
  use cands, only: ICSM, ISSM
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, LSM, n, nTU, nTUVX
  real(kind=wp), intent(in) :: Cin(*), TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout) :: Sout(*)
  integer(kind=iwp) :: IREFSM_Keep, nD
  real(kind=wp), allocatable :: CT(:), ST(:)

  IREFSM_Keep = IREFSM
  IREFSM = LSM
  ICSM = LSM
  ISSM = LSM
  nD = CI_Lib_Space(iState)%nDet(IREFSM)
  call mma_allocate(CT,max(n,nD),Label='CT')
  call mma_allocate(ST,max(n,nD),Label='ST')
  CT(1:n) = Cin(1:n)
  ST(:) = Zero
  call csdtvc(CT,ST,1,LSM,0)
  CT(1:nD) = ST(1:nD)
  call Lucia_Util('Sigma',CI_Vector=CT(:),Sigma_Vector=ST(:),nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)
  Sigma_on_disk = .true.
  call CSDTVC(Sout,ST,2,LSM,0)
  call mma_deallocate(ST)
  call mma_deallocate(CT)
  IREFSM = IREFSM_Keep
  ICSM = IREFSM
  ISSM = IREFSM

end subroutine Lucia_Sigma_Ref

!> A weighted checksum of A (the order matters)
function CheckSum(A,n) result(Chk)

  use Constants, only: Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: A(n)
  real(kind=wp) :: Chk
  integer(kind=iwp) :: i

  Chk = Zero
  do i=1,n
    Chk = Chk+abs(A(i))/real(i,kind=wp)
  end do

end function CheckSum

!> The CI vector CI (basis of the selected library, symmetry iSym of state
!> iState) for the active orbitals rotated by TU (per symmetry the square
!> nAsh x nAsh blocks): the transformation of the core (CI_Lib_Transform of
!> CI_Lib_Core). Under MOLCAS_CI_VERIFY=SIGMA compared with Lucia (SGUGA,
!> Faroald) and with SGUGA (the other libraries).
subroutine CI_Lib_Transform(iState,iSym,nTU,TU,n,CI)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, CI_File_To_Basis, CI_Lucia_To_Basis, CI_Faro_CI, CI_SG_Avail, &
                           CI_SG_CI, CI_Tim_Off, CI_Tim_On, CI_Tim_TraPsi, CI_Chk_TraPsi_Faro, CI_Chk_TraPsi
  use CI_Verify_Options, only: CI_Verify
  use Lucia_Interface, only: Lucia_Util
  use sguga, only: sg_traci
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, n
  real(kind=wp), intent(in) :: TU(nTU)
  real(kind=wp), intent(inout) :: CI(n)
  real(kind=wp) :: Dev, TimT(2)
  logical(kind=iwp) :: Ref_Lucia, Verify_SG
  real(kind=wp), allocatable :: CI_Ref(:), CI_SGUGA(:), ctemp(:)

  call CI_Tim_On(TimT)
  Ref_Lucia = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (CI_SG_CI() .or. CI_Faro_CI())
  Verify_SG = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. CI_SG_Avail() .and. (.not. CI_SG_CI())
  if (Ref_Lucia) then
    ! the reference of Lucia, in the Lucia basis
    call mma_allocate(CI_Ref,n,Label='CI_Lucia')
    call CI_Basis_To_Lucia(iState,iSym,n,CI,CI_Ref)
    call Lucia_Util('TRACI',CI_Vector=CI_Ref(:),TU=TU(:),nTU=nTU)
  end if
  if (Verify_SG) then
    ! the reference of SGUGA, in the order of the file
    call mma_allocate(CI_SGUGA,n,Label='CI_SGUGA')
    call CI_Basis_To_File(iState,iSym,n,CI,CI_SGUGA)
    call SG_TRACI(iState,iSym,TU,nTU,CI_SGUGA,n)
  end if

  call Raw_Transform(iState,iSym,nTU,TU,n,CI)

  if (Ref_Lucia) then
    call mma_allocate(ctemp,n,Label='ctemp')
    call CI_Lucia_To_Basis(iState,iSym,n,CI_Ref,ctemp)
    Dev = maxval(abs(CI(:)-ctemp(:)))
    if (CI_SG_CI()) then
      call CI_Chk_Report_Dev('CI_Lib_Transform (SGCI): SGUGA CI transformation vs Lucia',Dev,1.0e-10_wp,CI_Chk_TraPsi)
    else
      call CI_Chk_Report_Dev('CI_Lib_Transform (FARO): Faroald CI transformation vs Lucia',Dev,1.0e-10_wp,CI_Chk_TraPsi_Faro)
    end if
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_Ref)
  end if
  if (Verify_SG) then
    call mma_allocate(ctemp,n,Label='ctemp')
    call CI_File_To_Basis(iState,iSym,n,CI_SGUGA,ctemp)
    call CI_Chk_Report_Dev('CI_Lib_Transform: SGUGA CI transformation vs the CI library',maxval(abs(CI(:)-ctemp(:))), &
                    1.0e-10_wp*max(One,maxval(abs(CI(:)))),CI_Chk_TraPsi)
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_SGUGA)
  end if
  call CI_Tim_Off(CI_Tim_TraPsi,TimT)

end subroutine CI_Lib_Transform

!***********************************************************************************************************************************

!> The explicit block <I|H|J> + ECore of the core (CI_Lib_Block of
!> CI_Lib_Core, see there for the arguments). Under MOLCAS_CI_VERIFY=SIGMA,
!> once per run: the block in use is compared with the other libraries (with
!> FARO in the Faroald basis: columns H e_J of the sigma vectors), SGUGA's
!> own configuration list and selection are checked against it, and some of
!> its columns against the sigma vector (CI_Lib_Sigma); the CPU time of each
!> library is reported.
subroutine CI_Lib_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,nHP,HP)

  use CI_Lib_Checks, only: CI_Faro_CI, CI_Chk_IHJ, CI_SG_Avail, CI_SG_CI, CI_Tim_IHJ, CI_Tim_Off, CI_Tim_On
  use CI_Verify_Options, only: CI_Verify
  use CI_SG_Lucia_Bridge, only: CI_SG_Verify_Lucia_Walk
  use CI_Lib_Info, only: CI_Lib_Space
  use Definitions, only: wp, iwp

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nPMax, nHP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  integer(kind=iwp), intent(out) :: nP, IdxP(nPMax)
  real(kind=wp), intent(out) :: HP(nHP)
  real(kind=wp) :: T0, T1, TimT(2), W
  logical(kind=iwp) :: Verify

  call CI_Tim_On(TimT)
  Verify = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (.not. CI_Chk_IHJ)
  ! the walks of the Lucia CSFs, once per run
  if (Verify .and. (.not. CI_SG_CI()) .and. CI_SG_Avail() .and. (.not. CI_Faro_CI())) &
    call CI_SG_Verify_Lucia_Walk(iState,iSym,CI_Lib_Space(iState)%nCSF(iSym))
  call CWTime(T0,W)
  call Raw_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,nHP,HP)
  call CWTime(T1,W)
  if (Verify) then
    CI_Chk_IHJ = .true.
    call Check_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,HP,T1-T0)
  end if
  call CI_Tim_Off(CI_Tim_IHJ,TimT)

end subroutine CI_Lib_Block

!> CI development (SIGMA): the block HP of the library in use (CPU time
!> TBlk) against the other libraries, SGUGA's own selection
!> (Verify_SG_Explicit) and columns of the sigma vector (Block_Columns).
subroutine Check_Block(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,IdxP,HP,TBlk)

  use CI_Lib_Checks, only: CI_Faro_CI, CI_SG_Avail, CI_SG_CI
  use CI_Lib_Info, only: CI_Lib_Space
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nPMax, nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, HP(nTri_Elem(nP)), TBlk
  integer(kind=iwp) :: k, nCnf
  integer(kind=iwp), allocatable :: iCSFCnf(:), nCSFCnf(:), nSelCnf(:)
  real(kind=wp), allocatable :: HRef(:)
  logical(kind=iwp) :: SG_OK, UseSG
  character(len=8) :: Lab

  ! SGCI: the block of SGUGA, checked against columns of the sigma vector
  ! (itself compared with Lucia there)
  if (CI_SG_CI()) then
    call Block_Columns(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HP)
    return
  end if

  ! SGUGA is available without GAS (or through the relaxed space)
  SG_OK = CI_SG_Avail()
  ! SGEX: the block of the core is SGUGA's for RASSCF's Lucia
  UseSG = CI_Lib_SGEX(iState)
  if (UseSG) then
    Lab = 'SGUGA'
  else if (CI_Faro_CI()) then
    Lab = 'Faroald'
  else
    Lab = 'Lucia'
  end if
  write(u6,'(1X,3A,I0,A,F10.3,A)') 'CI_Lib_Block: ',trim(Lab),' block, ',nP,' CSFs, CPU ',TBlk,' s'

  ! the number of CSFs of each selected configuration
  nCnf = CI_Lib_nCnf(iState,iSym)
  call mma_allocate(nCSFCnf,max(1,nCnf),Label='nCSFCnf')
  call mma_allocate(iCSFCnf,max(1,CI_Lib_Space(iState)%nCSF(iSym)),Label='iCSFCnf')
  call CI_Lib_Cnf_Layout(iState,iSym,nCnf,nCSFCnf,size(iCSFCnf),iCSFCnf)
  call mma_allocate(nSelCnf,max(1,nPC),Label='nSelCnf')
  do k=1,nPC
    nSelCnf(k) = nCSFCnf(IPC(k))
  end do
  call mma_deallocate(iCSFCnf)
  call mma_deallocate(nCSFCnf)

  call mma_allocate(HRef,nTri_Elem(nP),Label='HRef')
  ! (with FARO the block is in the basis of Faroald: the references of
  ! Lucia and SGUGA are their columns H e_J in that basis, CI_Chk_Ref_Sigma)
  if (UseSG .and. (.not. CI_Faro_CI())) call Block_Ref('Lucia',Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP, &
    HP,HRef)
  if (CI_Faro_CI()) then
    call Block_Ref('Faroald H e_J',Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP,HP,HRef)
    call Block_Ref('Lucia H e_J',Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP,HP,HRef)
    if (SG_OK) call Block_Ref('SGUGA H e_J',Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP,HP,HRef)
  end if
  if (SG_OK) then
    if ((.not. UseSG) .and. (.not. CI_Faro_CI())) &
      call Block_Ref('SGUGA',Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP,HP,HRef)
    call Verify_SG_Explicit(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,HP)
  else
    write(u6,'(1X,A)') 'CI_Lib_Block: SGUGA block skipped (not available with GAS)'
  end if
  call mma_deallocate(HRef)
  call mma_deallocate(nSelCnf)

  call Block_Columns(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HP)

end subroutine Check_Block

!> CI development (SIGMA): the reference block HRef of LabRef, compared with
!> the block HP of the library in use (Lab).
subroutine Block_Ref(LabRef,Lab,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nSelCnf,nP,IdxP,HP,HRef)

  use Index_Functions, only: nTri_Elem
  use Definitions, only: wp, iwp, u6

  character(len=*), intent(in) :: LabRef, Lab
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nSelCnf(nPC), nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, HP(nTri_Elem(nP))
  real(kind=wp), intent(out) :: HRef(nTri_Elem(nP))
  real(kind=wp) :: T2, T3, W

  call CWTime(T2,W)
  select case (LabRef)
    case ('Lucia')
      ! (with SGEX: the Lucia block of the core)
      call CI_Lucia_Block(iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nP,HRef)
    case ('Faroald H e_J')
      call Block_Col_Ref('Faroald',iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HRef)
    case ('Lucia H e_J')
      call Block_Col_Ref('Lucia',iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HRef)
    case ('SGUGA H e_J')
      call Block_Col_Ref('SGUGA',iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HRef)
    case ('SGUGA')
      call CI_Lucia_Block_SG(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,nSelCnf,nP,IdxP,HRef)
  end select
  call CWTime(T3,W)
  write(u6,'(1X,3A,F10.3,A)') 'CI_Lib_Block: ',LabRef,' block, CPU ',T3-T2,' s'
  call Block_Compare(trim(Lab),LabRef,nP,IdxP,HP,HRef)

end subroutine Block_Ref

!> CI development (SIGMA): the reference block H(I,J) = (H e_J)(I) + ECore
!> delta_IJ from the sigma vectors of the unit vectors e_J, J over the CSFs
!> IdxP of the block: Faroald (the core sigma vector of the library in use,
!> with FARO), Lucia or SGUGA (CI_Chk_Ref_Sigma, in the basis of the library).
subroutine Block_Col_Ref(LabRef,iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,H)

  use CI_Lib_Checks, only: CI_Chk_Ref_Sigma
  use CI_Lib_Info, only: CI_Lib_Space
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp

  character(len=*), intent(in) :: LabRef
  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore
  real(kind=wp), intent(out) :: H(nTri_Elem(nP))
  integer(kind=iwp) :: I, J, n
  real(kind=wp), allocatable :: EJ(:), HEJ(:)

  n = CI_Lib_Space(iState)%nCSF(iSym)
  call mma_allocate(EJ,n,Label='EJ')
  call mma_allocate(HEJ,n,Label='HEJ')
  EJ(:) = Zero
  do J=1,nP
    EJ(IdxP(J)) = One
    if (LabRef == 'Faroald') then
      call Raw_Sigma(iState,0,iSym,iSym,.false.,.true.,nTU,TU,nTUVX,TUVX,n,EJ,n,HEJ)
    else
      call CI_Chk_Ref_Sigma(LabRef,iState,iSym,n,EJ,HEJ,nTU,TU,nTUVX,TUVX)
    end if
    EJ(IdxP(J)) = Zero
    do I=J,nP
      H(nTri_Elem(I-1)+J) = HEJ(IdxP(I))
    end do
    H(nTri_Elem(J)) = H(nTri_Elem(J))+ECore
  end do
  call mma_deallocate(HEJ)
  call mma_deallocate(EJ)

end subroutine Block_Col_Ref

!> CI development (SIGMA): the block HP (library Lab) against the reference
!> block HRef (LabRef), elementwise.
subroutine Block_Compare(Lab,LabRef,nP,IdxP,HP,HRef)

  use Index_Functions, only: nTri_Elem
  use Constants, only: Zero
  use Definitions, only: wp, iwp, u6

  character(len=*), intent(in) :: Lab, LabRef
  integer(kind=iwp), intent(in) :: nP, IdxP(nP)
  real(kind=wp), intent(in) :: HP(nTri_Elem(nP)), HRef(nTri_Elem(nP))
  real(kind=wp), parameter :: Thr = 1.0e-10_wp
  integer(kind=iwp) :: I, IJ, J
  real(kind=wp) :: ErrMax

  ErrMax = Zero
  do I=1,nP
    do J=1,I
      IJ = nTri_Elem(I-1)+J
      ErrMax = max(ErrMax,abs(HP(IJ)-HRef(IJ)))
      if (abs(HP(IJ)-HRef(IJ)) > Thr) then
        write(u6,*) 'CI_Lib_Block: <I|H|J> of ',Lab,' differs from ',LabRef
        write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IdxP(I),IdxP(J)
        write(u6,*) Lab,' =',HP(IJ)
        write(u6,*) LabRef,' =',HRef(IJ)
        call Abend()
      end if
    end do
  end do
  write(u6,'(1X,5A,I0,A,ES10.3)') 'CI_Lib_Block: ',Lab,' verified against ',LabRef,', ',nP,' CSFs, max. error ',ErrMax

end subroutine Block_Compare

!> CI development (SIGMA): the block HP against columns H e_J of the sigma
!> vector (CI_Lib_Sigma, checked itself) for at most mxCol CSFs J of the
!> block, evenly spaced.
subroutine Block_Columns(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nP,IdxP,HP)

  use Lucia_data, only: Sigma_on_disk
  use CI_Lib_Info, only: CI_Lib_Space
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nP, IdxP(nP)
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, HP(nTri_Elem(nP))
  integer(kind=iwp), parameter :: mxCol = 32
  real(kind=wp), parameter :: Thr = 1.0e-10_wp
  integer(kind=iwp) :: I, J, JJ, n, nCol
  real(kind=wp) :: Err, ErrMax, Ref
  logical(kind=iwp) :: Sigma_Save
  real(kind=wp), allocatable :: EJ(:), HEJ(:)

  n = CI_Lib_Space(iState)%nCSF(iSym)
  call mma_allocate(EJ,n,Label='EJ')
  call mma_allocate(HEJ,n,Label='HEJ')
  Sigma_Save = Sigma_on_disk
  nCol = min(nP,mxCol)
  ErrMax = Zero
  do JJ=1,nCol
    J = 1+((JJ-1)*(nP-1))/max(nCol-1,1)
    EJ(:) = Zero
    EJ(IdxP(J)) = One
    call CI_Lib_Sigma(iState,0,iSym,iSym,.false.,.true.,nTU,TU,nTUVX,TUVX,n,EJ,n,HEJ)
    do I=1,nP
      Ref = HEJ(IdxP(I))
      if (I == J) Ref = Ref+ECore
      Err = abs(HP(nTri_Elem(max(I,J)-1)+min(I,J))-Ref)
      ErrMax = max(ErrMax,Err)
      if (Err > Thr) then
        write(u6,*) 'CI_Lib_Block: <I|H|J> differs from (H e_J)(I)'
        write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IdxP(I),IdxP(J)
        write(u6,*) 'CI_Lib_Block   =',HP(nTri_Elem(max(I,J)-1)+min(I,J))
        write(u6,*) 'CI_Lib_Sigma =',Ref
        call Abend()
      end if
    end do
  end do
  Sigma_on_disk = Sigma_Save
  write(u6,'(1X,A,I0,A,I0,A,ES10.3)') 'CI_Lib_Block verified against CI_Lib_Sigma: ',nCol,' of ',nP,' columns, max. error ',ErrMax
  call mma_deallocate(HEJ)
  call mma_deallocate(EJ)

end subroutine Block_Columns

!> CI development (SIGMA): SGUGA's own configuration list and explicit
!> subspace (SG_Mk_Cnf_List, SG_Select_Explicit, SG_Explicit_Block) against
!> the library in use: (1) each configuration of the library is exactly one
!> SGUGA configuration (same CSFs), (2) the two selections are compared (the
!> keys differ: the configuration average of the library, the mean exact CSF
!> diagonal of SGUGA), (3) the SGUGA block over the configurations IPC of the
!> library has the eigenvalues of the block in use, HP (nP CSFs).
subroutine Verify_SG_Explicit(iState,iSym,nTU,TU,nTUVX,TUVX,ECore,nPC,IPC,nPMax,nP,HP)

  use CI_Lib_Checks, only: CI_Basis_To_File
  use CI_Lib_Info, only: CI_Lib_Space
  use sguga, only: SGCnfList, SG_Cnf_Avg_Diag, SG_Explicit_Block, SG_Free_Cnf_List, SG_Mk_Cnf_List, SG_Select_Explicit
  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  use Definitions, only: wp, iwp, u6

  integer(kind=iwp), intent(in) :: iState, iSym, nTU, nTUVX, nPC, IPC(nPC), nPMax, nP
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX), ECore, HP(nTri_Elem(nP))
  real(kind=wp), parameter :: Thr = 1.0e-10_wp
  type(SGCnfList) :: Cnf
  integer(kind=iwp) :: I, iOffL, k, kG, nBoth, nConf, nGCSF, nGSel, nLSel
  integer(kind=iwp), allocatable :: CnfOfG(:), GCSF(:), GSel(:), GtoL(:), LtoG(:), LtoGCnf(:), nCSFL(:)
  logical(kind=iwp), allocatable :: SelG(:), SelL(:)
  real(kind=wp) :: Dev
  real(kind=wp), allocatable :: DSG(:), EG(:), EL(:), HG(:), Idx(:), Map(:)

  nConf = CI_Lib_Space(iState)%nCSF(iSym)
  call SG_Mk_Cnf_List(iState,iSym,Cnf)
  ! library CSF -> split-graph CSF, and the SGUGA configuration of each
  call mma_allocate(Idx,nConf,Label='Idx')
  call mma_allocate(Map,nConf,Label='Map')
  do I=1,nConf
    Idx(I) = real(I,kind=wp)
  end do
  call CI_Basis_To_File(iState,iSym,nConf,Idx,Map)
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
  ! (1) configurations of the library -> SGUGA configurations, a bijection
  call mma_allocate(nCSFL,CI_Lib_nCnf(iState,iSym),Label='nCSFL')
  call mma_allocate(GtoL,nConf,Label='GtoL')
  call CI_Lib_Cnf_Layout(iState,iSym,size(nCSFL),nCSFL,nConf,GtoL)
  call mma_deallocate(GtoL)
  if (size(nCSFL) /= Cnf%nCnf) then
    write(u6,*) 'CI_Lib_Block: number of configurations, CI back end and SGUGA:',size(nCSFL),Cnf%nCnf
    call Abend()
  end if
  call mma_allocate(LtoGCnf,size(nCSFL),Label='LtoGCnf')
  call mma_allocate(SelG,Cnf%nCnf,Label='SelG')
  SelG(:) = .false.
  iOffL = 0
  do k=1,size(nCSFL)
    kG = CnfOfG(LtoG(iOffL+1))
    if (any(CnfOfG(LtoG(iOffL+1:iOffL+nCSFL(k))) /= kG) .or. (nCSFL(k) /= Cnf%iOff(kG+1)-Cnf%iOff(kG)) .or. SelG(kG)) then
      write(u6,*) 'CI_Lib_Block: configuration',k,' of the CI back end is not one SGUGA configuration'
      call Abend()
    end if
    SelG(kG) = .true.
    LtoGCnf(k) = kG
    iOffL = iOffL+nCSFL(k)
  end do
  write(u6,'(1X,A,I0,A)') 'CI_Lib_Block: SGUGA configuration list verified, ',Cnf%nCnf,' configurations as in the CI back end'
  ! (2) the SGUGA selection with its configuration-average diagonal (the key
  ! of SGCI and of Lucia: differences come from ties only)
  call mma_allocate(DSG,nConf,Label='DSG')
  call SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,DSG)
  call CI_HEXS_Proj(CI_HEXS_Keep_SG,DSG)
  call mma_allocate(GSel,Cnf%nCnf,Label='GSel')
  call mma_allocate(GCSF,max(1,nPMax),Label='GCSF')
  call SG_Select_Explicit(Cnf,DSG,nPMax,nGSel,GSel,nGCSF,GCSF)
  call mma_deallocate(DSG)
  call mma_allocate(SelL,Cnf%nCnf,Label='SelL')
  SelG(:) = .false.
  SelL(:) = .false.
  SelG(GSel(1:nGSel)) = .true.
  SelL(LtoGCnf(IPC(1:nPC))) = .true.
  nBoth = count(SelG .and. SelL)
  nLSel = count(SelL)
  write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'CI_Lib_Block: SGUGA selection (configuration average): ',nGSel,' configurations, ', &
                                              nGCSF,' CSFs; CI back end: ',nLSel,' configurations; ',nBoth,' in common, ', &
                                              count(SelG .neqv. SelL),' differ'
  ! (3) SGUGA block over the configurations selected by the library (for GAS
  ! the matrix elements in the relaxed space)
  call mma_allocate(GtoL,max(1,nPC),Label='GtoL')
  GtoL(1:nPC) = LtoGCnf(IPC(1:nPC))
  call mma_allocate(HG,nTri_Elem(nP),Label='HG')
  if (CI_Lib_Space(iState)%iRel > 0) then
    call CI_SG_Rel_Map(iState,Cnf%LSM)
    call SG_Explicit_Block(iState,Cnf,nPC,GtoL,ECore,TU,nTU,TUVX,nTUVX,HG,iRel=CI_Lib_Space(iState)%iRel, &
      RelMap=CI_Rel_Map(Cnf%LSM)%M)
  else
    call SG_Explicit_Block(iState,Cnf,nPC,GtoL,ECore,TU,nTU,TUVX,nTUVX,HG)
  end if
  call mma_allocate(EG,nP,Label='EG')
  call mma_allocate(EL,nP,Label='EL')
  call Block_Eigenvalues(HG,nP,EG)
  HG(:) = HP(:)
  call Block_Eigenvalues(HG,nP,EL)
  Dev = maxval(abs(EG(:)-EL(:)))
  write(u6,'(1X,A,I0,A,ES10.3)') 'CI_Lib_Block: SG_Explicit_Block, eigenvalues verified against the block in use, ',nP, &
                                 ' CSFs, max. error ',Dev
  if (Dev > Thr*max(One,maxval(abs(EL(:))))) then
    write(u6,*) 'CI_Lib_Block: eigenvalues of the SGUGA explicit block differ'
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

!> Eigenvalues E (increasing) of the lower-triangle packed symmetric matrix
!> A, which is overwritten.
subroutine Block_Eigenvalues(A,n,E)

  use Index_Functions, only: nTri_Elem
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Zero
  use Definitions, only: wp, iwp

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

end module CI_Lib
