!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.             *
!                                                                      *
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module UMP2_Global

use Definitions, only: wp, iwp
use Constants, only: Zero

implicit none
private

! ---------------------------------------------------------------------
! Shared state for the open-shell MP2 energy implementations.
!
! The C1 UMP2 driver supports conventional stored AO integrals and
! Cholesky vectors in serial and MPI calculations.  Non-C1 references
! are handled by the symmetry-blocked open-shell MP2 driver.
!
! Para_Info remains the source of MPI rank information.  The C1
! Cholesky driver assembles complete replicated transformed vectors
! before the occupied-pair contractions are distributed across ranks.
! ---------------------------------------------------------------------

! Dimensions of the reference.
!
! nBas: number of AO basis functions.
! nOrb: number of retained MOs per spin, excluding basis directions
!       removed during SCF.
!
! These are scalars for the C1 drivers.  The symmetry-blocked path
! maintains its irrep-resolved dimensions separately.

integer(kind=iwp) :: nSym = 0
integer(kind=iwp) :: nBas = 0
integer(kind=iwp) :: nOrb = 0

! Total occupied-orbital counts, including frozen core.

integer(kind=iwp) :: nOccA = 0
integer(kind=iwp) :: nOccB = 0

! Number of occupied orbitals frozen in each spin space.
!
! Required:
!   0 <= nFro <= min(nOccA,nOccB)
!
! ExplicitFrozen distinguishes an explicit "Frozen = 0" request
! from an input in which no frozen-core choice was supplied.

integer(kind=iwp) :: nFro = 0
logical(kind=iwp) :: ExplicitFrozen = .false.

! Correlated occupied and virtual dimensions.
!
! nOA = nOccA - nFro
! nOB = nOccB - nFro
! nVA = nOrb  - nOccA
! nVB = nOrb  - nOccB

integer(kind=iwp) :: nOA = 0
integer(kind=iwp) :: nOB = 0
integer(kind=iwp) :: nVA = 0
integer(kind=iwp) :: nVB = 0

! Input/reference information.
!
! These flags must be populated and checked by the driver.
! Their default values do not certify a valid calculation.

integer(kind=iwp) :: iPL = 0
character(len=8) :: ReferenceMethod = '        '
logical(kind=iwp) :: DoCholesky = .false.

! Minimum allowed absolute orbital-energy denominator, in hartree.
!
! This is an error-detection threshold, not a denominator shift.
! The driver may pass it directly to UMP2_Energy.

real(kind=wp) :: DenTol = 1.0e-12_wp

! Full canonical MO coefficients.
!
! Allocate:
!   CAlpha(nBas,nOrb)
!   CBeta (nBas,nOrb)
!
! First index: AO basis function.
! Second index: canonical MO.
!
! Occupied orbitals precede virtual orbitals in each spin set.
! The reference reader must preserve the correspondence between
! coefficient columns and orbital energies.

real(kind=wp), allocatable :: CAlpha(:,:)
real(kind=wp), allocatable :: CBeta(:,:)

! Full canonical orbital energies, including frozen occupied orbitals.
!
! Allocate:
!   EOrbA(nOrb)
!   EOrbB(nOrb)

real(kind=wp), allocatable :: EOrbA(:)
real(kind=wp), allocatable :: EOrbB(:)

! Correlated occupied and virtual orbital energies.
!
! Allocate:
!   EOccA(nOA), EVirA(nVA)
!   EOccB(nOB), EVirB(nVB)
!
! Space construction:
!   EOccA = EOrbA(nFro+1:nOccA)
!   EOccB = EOrbB(nFro+1:nOccB)
!   EVirA = EOrbA(nOccA+1:nOrb)
!   EVirB = EOrbB(nOccB+1:nOrb)
!
! Zero-sized occupied or virtual spaces are permitted.

real(kind=wp), allocatable :: EOccA(:)
real(kind=wp), allocatable :: EOccB(:)
real(kind=wp), allocatable :: EVirA(:)
real(kind=wp), allocatable :: EVirB(:)

! Energies in hartree.
!
! EAA, EBB, EAB are signed correlation contributions:
!   ECorr = EAA + EBB + EAB
!   ETotal = ESCF + ECorr
!
! EnergyReady must remain false unless evaluation succeeds.
! Zero-initialized energies must not be treated as valid results.

real(kind=wp) :: ESCF = Zero
real(kind=wp) :: EAA = Zero
real(kind=wp) :: EBB = Zero
real(kind=wp) :: EAB = Zero
real(kind=wp) :: ECorr = Zero
real(kind=wp) :: ETotal = Zero

logical(kind=iwp) :: EnergyReady = .false.

public :: nSym, nBas, nOrb, nOccA, nOccB, nFro
public :: nOA, nOB, nVA, nVB
public :: iPL, ReferenceMethod, DoCholesky, ExplicitFrozen
public :: DenTol
public :: CAlpha, CBeta, EOrbA, EOrbB
public :: EOccA, EOccB, EVirA, EVirB
public :: ESCF, EAA, EBB, EAB, ECorr, ETotal, EnergyReady
public :: UMP2_Clean

contains

subroutine UMP2_Clean()

  use stdalloc, only: mma_deallocate

  ! All module arrays must be allocated through mma_allocate.
  ! The safe option permits cleanup before all arrays exist.

  call mma_deallocate(CAlpha,safe='*')
  call mma_deallocate(CBeta,safe='*')

  call mma_deallocate(EOrbA,safe='*')
  call mma_deallocate(EOrbB,safe='*')

  call mma_deallocate(EOccA,safe='*')
  call mma_deallocate(EOccB,safe='*')
  call mma_deallocate(EVirA,safe='*')
  call mma_deallocate(EVirB,safe='*')

  ! Reset state so a subsequent invocation cannot reuse old values.

  nSym = 0
  nBas = 0
  nOrb = 0

  nOccA = 0
  nOccB = 0
  nFro = 0

  nOA = 0
  nOB = 0
  nVA = 0
  nVB = 0

  iPL = 0
  ReferenceMethod = '        '
  DoCholesky = .false.
  ExplicitFrozen = .false.

  DenTol = 1.0e-12_wp

  ESCF = Zero
  EAA = Zero
  EBB = Zero
  EAB = Zero
  ECorr = Zero
  ETotal = Zero

  EnergyReady = .false.

end subroutine UMP2_Clean

end module UMP2_Global