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
!> generic CI operations CI_Lib_* of the core (CI_Operations, see there),
!> with the checks of MOLCAS_CI_VERIFY (CI_Lib_Checks) around those that
!> are verified. The checks apply to the CI state set up by
!> CI_Lib_Initialize from the CI input (CI_Lib_Input_State); a state the
!> program set up itself gets the core operations unchanged.
module CI_Library

use CI_Operations, Raw_Dens => CI_Lib_Dens, Raw_Transform => CI_Lib_Transform

implicit none
private :: Check_T1DM, Raw_Dens, Raw_Transform

public :: CI_Lib_Dens, CI_Lib_Transform

contains

!> The 1- and 2-particle (transition) densities of the core (CI_Lib_Dens of
!> CI_Operations, see there for the layouts). Under MOLCAS_CI_VERIFY=SIGMA a
!> one-body transition density of one symmetry (spin-free, not packed) is
!> compared, together with its spin part, with Lucia and SGUGA (those not
!> in use).
subroutine CI_Lib_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Packed)

  use CI_Lib_Checks, only: Tim_Off, Tim_On, Tim_PDM, Tim_T1DM
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
  call Tim_On(TimT)
  call Raw_Dens(iState,I12,iOpSpin,iLSym,iRSym,nL,L,nR,R,n1,Rho1,n2,Rho2,Packed)
  if (CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (.not. Pk) .and. (I12 == 1) .and. (iOpSpin == 0) .and. &
      (iLSym == iRSym) .and. (nL == nR)) call Check_T1DM(iState,iLSym,nL,L,R,n1,Rho1,n2)
  if (Pk) then
    call Tim_Off(Tim_PDM,TimT)
  else
    call Tim_Off(Tim_T1DM,TimT)
  end if

end subroutine CI_Lib_Dens

!> CI development (SIGMA): the one-body transition density Rho1 (spin-free,
!> layout of CI_Lib_Dens) of Bra and Ket and its spin part (computed here)
!> against Lucia (no spin density for singlets) and SGUGA (sg_trans_pdm,
!> whose layout is the transpose), those not in use; in the level layout
!> T(p+(q-1)*nAct) = <Bra|E_qp|Ket> (Lucia's).
subroutine Check_T1DM(iState,LSM,nVec,Bra_Vec,Ket_Vec,n1,Rho1,n2)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, SG_Avail, SG_Trans_PDM_G, T1DM_Faro_Verified, T1DM_SG_Verified
  use CI_Info, only: CI_St
  use CI_Op_Defs, only: CI_Lib, CI_Lib_Lucia, CI_Lib_SGUGA
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

  NAC = CI_St(iState)%nLev
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
  if (CI_Lib /= CI_Lib_Lucia) then
    call CI_Basis_To_Lucia(iState,LSM,nVec,Bra_Vec,BL)
    call CI_Basis_To_Lucia(iState,LSM,nVec,Ket_Vec,KL)
    call mma_allocate(DTMP,NAC*NAC,Label='DTMP')
    call mma_allocate(DSTMP,NAC*NAC,Label='DSTMP')
    call Lucia_Util('Densi',CI_Vector=KL(:),RVec=BL(:))
    Dev(1) = maxval(abs(DTMP(1:NAC*NAC)-T1DM(1:NAC*NAC)))
    Dev(2) = Zero
    if (CI_St(iState)%iSpin > 1) Dev(2) = maxval(abs(DSTMP(1:NAC*NAC)-T1SDM(1:NAC*NAC)))
    Done = T1DM_Faro_Verified
    call Report_Dev('Mk_T1DM: T1DM, T1SDM vs Lucia',maxval(Dev),ThrD,Done)
    call mma_deallocate(DSTMP)
    call mma_deallocate(DTMP)
  end if
  if ((CI_Lib /= CI_Lib_SGUGA) .and. SG_Avail()) then
    call CI_Basis_To_File(iState,LSM,nVec,Bra_Vec,BL)
    call CI_Basis_To_File(iState,LSM,nVec,Ket_Vec,KL)
    call mma_allocate(TS,NAC*NAC,Label='TS')
    call mma_allocate(TSS,NAC*NAC,Label='TSS')
    call SG_Trans_PDM_G(iState,BL,KL,nVec,LSM,NAC,TS,TSS,.true.)
    Dev(:) = Zero
    do q=1,NAC
      do p=1,NAC
        Dev(1) = max(Dev(1),abs(T1DM(p+(q-1)*NAC)-TS(q+(p-1)*NAC)))
        if (CI_Lib /= CI_Lib_Lucia) Dev(2) = max(Dev(2),abs(T1SDM(p+(q-1)*NAC)-TSS(q+(p-1)*NAC)))
      end do
    end do
    Done = T1DM_SG_Verified
    call Report_Dev('Mk_T1DM: T1DM, T1SDM vs SGUGA',maxval(Dev),ThrD,Done)
    call mma_deallocate(TSS)
    call mma_deallocate(TS)
  end if
  T1DM_Faro_Verified = .true.
  T1DM_SG_Verified = .true.
  call mma_deallocate(KL)
  call mma_deallocate(BL)
  call mma_deallocate(T1SDM)
  call mma_deallocate(T1DM)

end subroutine Check_T1DM

!> The CI vector CI (basis of the selected library, symmetry iSym of state
!> iState) for the active orbitals rotated by TU (per symmetry the square
!> nAsh x nAsh blocks): the transformation of the core (CI_Lib_Transform of
!> CI_Operations). Under MOLCAS_CI_VERIFY=SIGMA compared with Lucia (SGUGA,
!> Faroald) and with SGUGA (the other libraries).
subroutine CI_Lib_Transform(iState,iSym,nTU,TU,n,CI)

  use CI_Lib_Checks, only: CI_Basis_To_File, CI_Basis_To_Lucia, CI_File_To_Basis, CI_Lucia_To_Basis, Faro_CI, SG_Avail, &
                           SG_CI, Tim_Off, Tim_On, Tim_TraPsi, TraPsi_Faro_Verified, TraPsi_Verified
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

  call Tim_On(TimT)
  Ref_Lucia = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. (SG_CI() .or. Faro_CI())
  Verify_SG = CI_Verify('SIGMA') .and. CI_Lib_Input_State(iState) .and. SG_Avail() .and. (.not. SG_CI())
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
    if (SG_CI()) then
      call Report_Dev('Mk_Tra_Psi (SGCI): SGUGA CI transformation vs Lucia',Dev,1.0e-10_wp,TraPsi_Verified)
    else
      call Report_Dev('Mk_Tra_Psi (FARO): Faroald CI transformation vs Lucia',Dev,1.0e-10_wp,TraPsi_Faro_Verified)
    end if
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_Ref)
  end if
  if (Verify_SG) then
    call mma_allocate(ctemp,n,Label='ctemp')
    call CI_File_To_Basis(iState,iSym,n,CI_SGUGA,ctemp)
    call Report_Dev('Mk_Tra_Psi: SGUGA CI transformation vs the CI library',maxval(abs(CI(:)-ctemp(:))), &
                    1.0e-10_wp*max(One,maxval(abs(CI(:)))),TraPsi_Verified)
    call mma_deallocate(ctemp)
    call mma_deallocate(CI_SGUGA)
  end if
  call Tim_Off(Tim_TraPsi,TimT)

end subroutine CI_Lib_Transform

end module CI_Library
