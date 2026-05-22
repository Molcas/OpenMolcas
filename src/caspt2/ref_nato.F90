!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine REF_NATO(DREF,nDREF,CMO,nCMO,OCC,nOcc,CNAT,nCNAT)
! Purpose: compute natural orbitals and natural occupation numbers
! for the reference wave function.
! Given DREF, a triangular density matrix
! in active/active MO basis, and symmetry-blocked array CMO of MO
! coefficients, return array of  natural occupation numbers and MO
! coefficients of  natural orbitals. Frozen, inactive and virtual
! orbitals are copied unchanged.

use caspt2_module, only: NASH, NBAS, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDREF, nCMO, nOcc, nCNAT
real(kind=wp), intent(in) :: DREF(nDREF), CMO(nCMO)
real(kind=wp), intent(out) :: OCC(nOcc), CNAT(nCNAT)
integer(kind=iwp) :: I, ICMO, IDREF, II, IOCC, ISYM, LIJ, NA, NB, NF, NFI, NI, NSD, NTMP
real(kind=wp) :: OC
real(kind=wp), allocatable :: TMP(:)

IDREF = 0
IOCC = 0
ICMO = 0
do ISYM=1,NSYM
  NF = NFRO(ISYM)
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NB = NBAS(ISYM)
  ! Frozen and inactive orbitals:
  NFI = NF+NI
  if (NFI > 0) then
    OCC(IOCC+1:IOCC+NFI) = Two
    IOCC = IOCC+NFI
    CNAT(ICMO+1:ICMO+NB*NFI) = CMO(ICMO+1:ICMO+NB*NFI)
    ICMO = ICMO+NB*NFI
  end if
  ! Active orbitals:
  if (NA > 0) then
    NTMP = (NA*(NA+1))/2
    call mma_allocate(TMP,NTMP,LABEL='TMP')
    CNAT(ICMO+1:ICMO+NB*NA) = CMO(ICMO+1:ICMO+NB*NA)
    ! For correct ordering, change sign.
    LIJ = 0
    do I=1,NA
      II = I+IDREF
      TMP(LIJ+1:LIJ+I) = -DREF((II*(II-1))/2+IDREF+1:(II*(II+1))/2)
      LIJ = LIJ+I
    end do
    call NIDiag(TMP,CNAT(ICMO+1),NA,NB)
    call JACORD(TMP,CNAT(ICMO+1),NA,NB)
    call VEIG(NA,TMP,OCC(IOCC+1))
    call mma_deallocate(TMP)
    ! Change back to positive sign.
    OCC(IOCC+1:IOCC+NA) = -OCC(IOCC+1:IOCC+NA)
    ! Certain CAS or RAS wave functions can legitimately have
    ! occupation numbers that are exactly 0 or 2. These may become
    ! inappropriate by rounding. Fix that as well.
    do I=1,NA
      OC = OCC(IOCC+I)
      if (OC < Zero) OC = Zero
      if (OC > Two) OC = Two
      OCC(IOCC+I) = OC
    end do
    IDREF = IDREF+NA
    IOCC = IOCC+NA
    ICMO = ICMO+NB*NA
  end if
  ! Secondary and deleted orbitals:
  NSD = NB-(NFI+NA)
  if (NSD > 0) then
    OCC(IOCC+1:IOCC+NSD) = Zero
    IOCC = IOCC+NSD
    CNAT(ICMO+1:ICMO+NB*NSD) = CMO(ICMO+1:ICMO+NB*NSD)
    ICMO = ICMO+NB*NSD
  end if
end do

end subroutine REF_NATO
