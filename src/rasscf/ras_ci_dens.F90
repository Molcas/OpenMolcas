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

!> RASSCF's densities from the CI library (CI_Library, CI_Lib_Dens), in
!> the layouts of RASSCF.
module RAS_CI_Dens

use CI_Library, only: CI_Lib_Dens, CI_Lib_Level_Map
use CI_Info, only: CI_St
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

public :: RAS_PDMs, RAS_T1DM

contains

!> The densities of the state CIVec (basis of the selected library) in the
!> layout of RASSCF: D, SD (triangular), P, PA (packed pairs), in level
!> order, from CI_Lib_Dens in its packed layout (any library).
subroutine RAS_PDMs(CIVec,nCIVEC,D,SD,P,PA,nD,nP)

  integer(kind=iwp), intent(in) :: nCIVEC, nD, nP
  real(kind=wp), intent(in) :: CIVEC(nCIVEC)
  real(kind=wp), intent(out), optional :: D(:), SD(:), P(:), PA(:)
  ! the CI state of RASSCF
  integer(kind=iwp), parameter :: iState = 1
  integer(kind=iwp) :: I12, LSM, NACPAR, NACPR2
  real(kind=wp), allocatable :: D_loc(:), R2(:), SD_loc(:)

  LSM = CI_St(iState)%STSym
  ! the packed active pairs and pairs of pairs
  NACPAR = nTri_Elem(CI_St(iState)%nLev)
  NACPR2 = nTri_Elem(NACPAR)
  I12 = 1
  if (present(P) .or. present(PA)) I12 = 2
  call mma_allocate(D_loc,max(1,NACPAR),Label='D_loc')
  call mma_allocate(R2,max(1,2*NACPR2),Label='R2')
  R2(:) = Zero
  call CI_Lib_Dens(iState,I12,0,LSM,LSM,nCIVEC,CIVEC,nCIVEC,CIVEC,size(D_loc),D_loc,size(R2),R2,Packed=.true.)
  if (present(D)) then
    D(1:nD) = Zero
    D(1:NACPAR) = D_loc(1:NACPAR)
  end if
  if (present(P)) P(1:nP) = R2(1:nP)
  if (present(PA)) PA(1:nP) = R2(NACPR2+1:NACPR2+nP)
  if (present(SD)) then
    call mma_allocate(SD_loc,max(1,NACPAR),Label='SD_loc')
    call CI_Lib_Dens(iState,1,1,LSM,LSM,nCIVEC,CIVEC,nCIVEC,CIVEC,size(SD_loc),SD_loc,size(R2),R2,Packed=.true.)
    SD(1:nD) = Zero
    SD(1:NACPAR) = SD_loc(1:NACPAR)
    call mma_deallocate(SD_loc)
  end if
  call mma_deallocate(R2)
  call mma_deallocate(D_loc)

end subroutine RAS_PDMs

!> The one-body transition density (and its spin part) of the CI vectors
!> Bra and Ket (basis of the selected library, the symmetry of the state)
!> in the layout of RASSCF: T1DM(p+(q-1)*NAC) = <Bra|E_qp|Ket> over the
!> levels (Lucia's layout), a permutation of Rho1 of CI_Lib_Dens.
subroutine RAS_T1DM(Bra_Vec,Ket_Vec,nVec,T1DM,nT1DM,T1SDM)

  integer(kind=iwp), intent(in) :: nVec, nT1DM
  real(kind=wp), intent(in) :: Bra_Vec(nVec), Ket_Vec(nVec)
  real(kind=wp), intent(out) :: T1DM(nT1DM)
  real(kind=wp), optional, intent(out) :: T1SDM(nT1DM)
  ! the CI state of RASSCF
  integer(kind=iwp), parameter :: iState = 1
  integer(kind=iwp) :: LSM, NAC
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: R1(:), R2(:)

  LSM = CI_St(iState)%STSym
  NAC = CI_St(iState)%nLev
  call mma_allocate(L2O,max(1,NAC),Label='L2O')
  call CI_Lib_Level_Map(iState,NAC,L2O)
  call mma_allocate(R1,max(1,NAC**2),Label='R1')
  call mma_allocate(R2,1,Label='R2')
  call CI_Lib_Dens(iState,1,0,LSM,LSM,nVec,Bra_Vec,nVec,Ket_Vec,size(R1),R1,size(R2),R2)
  call Level_Layout(NAC,L2O,R1,nT1DM,T1DM)
  if (present(T1SDM)) then
    call CI_Lib_Dens(iState,1,1,LSM,LSM,nVec,Bra_Vec,nVec,Ket_Vec,size(R1),R1,size(R2),R2)
    call Level_Layout(NAC,L2O,R1,nT1DM,T1SDM)
  end if
  call mma_deallocate(R2)
  call mma_deallocate(R1)
  call mma_deallocate(L2O)

end subroutine RAS_T1DM

!> T(p+(q-1)*NAC) = R1(L2O(q)+(L2O(p)-1)*NAC) over the levels p, q (the
!> rest of T zero)
subroutine Level_Layout(NAC,L2O,R1,nT,T)

  integer(kind=iwp), intent(in) :: NAC, L2O(NAC), nT
  real(kind=wp), intent(in) :: R1(NAC**2)
  real(kind=wp), intent(out) :: T(nT)
  integer(kind=iwp) :: p, q

  T(:) = Zero
  do q=1,NAC
    do p=1,NAC
      T(p+(q-1)*NAC) = R1(L2O(q)+(L2O(p)-1)*NAC)
    end do
  end do

end subroutine Level_Layout

end module RAS_CI_Dens
