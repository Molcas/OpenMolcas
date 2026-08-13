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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine NATORB_CASPT2(DMAT,nDMAT,CMO,nCMO,OCC,nOcc,CNAT,nCNAT)
! Given DMAT, symmetry-blocked array of triangular
! density matrices in MO basis, and symmetry-blocked
! array CMO of MO coefficients, return array of
! natural occupation numbers and MO coefficients of
! natural orbitals.

use Index_Functions, only: nTri_Elem
use caspt2_module, only: NBAS, NDEL, NFRO, NORB, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDMAT, nCMO, nOcc, nCNAT
real(kind=wp), intent(in) :: DMAT(nDMAT), CMO(nCMO)
real(kind=wp), intent(out) :: OCC(nOcc), CNAT(nCNAT)
integer(kind=iwp) :: ICMO, IDMAT, IOCC, ISYM, NB, ND, NF, NO, NTMP
real(kind=wp), allocatable :: TMP(:)

IDMAT = 0
IOCC = 0
ICMO = 0
do ISYM=1,NSYM
  NF = NFRO(ISYM)
  NO = NORB(ISYM)
  ND = NDEL(ISYM)
  NB = NBAS(ISYM)
  !  Frozen orbitals:
  if (NF > 0) then
    OCC(IOCC+1:IOCC+NF) = Two
    IOCC = IOCC+NF
    CNAT(ICMO+1:ICMO+NB*NF) = CMO(ICMO+1:ICMO+NB*NF)
    ICMO = ICMO+NB*NF
  end if
  ! Inactive, active, and secondary orbitals:
  if (NO > 0) then
    NTMP = nTri_Elem(NO)
    call mma_allocate(TMP,NTMP,Label='TMP')
    CNAT(ICMO+1:ICMO+NB*NO) = CMO(ICMO+1:ICMO+NB*NO)
    ! For correct order, change sign.
    TMP(:) = -DMAT(IDMAT+1:IDMAT+NTMP)
    call NIDiag(TMP,CNAT(ICMO+1),NO,NB)
    call JACORD(TMP,CNAT(ICMO+1),NO,NB)
    call VEIG(NO,TMP,OCC(IOCC+1))
    ! Change back to positive sign.
    OCC(IOCC+1:IOCC+NO) = -OCC(IOCC+1:IOCC+NO)
    IDMAT = IDMAT+NTMP
    IOCC = IOCC+NO
    ICMO = ICMO+NB*NO
    call mma_deallocate(TMP)
  end if
  ! Deleted orbitals:
  if (ND > 0) then
    OCC(IOCC+1:IOCC+ND) = Zero
    IOCC = IOCC+ND
    CNAT(ICMO+1:ICMO+NB*ND) = CMO(ICMO+1:ICMO+NB*ND)
    ICMO = ICMO+NB*ND
  end if
end do

end subroutine NATORB_CASPT2
