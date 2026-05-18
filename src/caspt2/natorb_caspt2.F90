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

use caspt2_module, only: NBAS, NDEL, NFRO, NORB, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
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
    call DCOPY_(NF,[Two],0,OCC(IOCC+1),1)
    IOCC = IOCC+NF
    call DCOPY_(NB*NF,CMO(ICMO+1),1,CNAT(ICMO+1),1)
    ICMO = ICMO+NB*NF
  end if
  ! Inactive, active, and secondary orbitals:
  if (NO > 0) then
    NTMP = (NO*(NO+1))/2
    call mma_allocate(TMP,NTMP,Label='TMP')
    call DCOPY_(NB*NO,CMO(ICMO+1),1,CNAT(ICMO+1),1)
    ! For correct order, change sign.
    call DYAX(NTMP,-One,DMAT(IDMAT+1),1,TMP,1)
    call NIDiag(TMP,CNAT(ICMO+1),NO,NB)
    call JACORD(TMP,CNAT(ICMO+1),NO,NB)
    call VEIG(NO,TMP,OCC(IOCC+1))
    ! Change back to positive sign.
    call DSCAL_(NO,-One,OCC(IOCC+1),1)
    IDMAT = IDMAT+NTMP
    IOCC = IOCC+NO
    ICMO = ICMO+NB*NO
    call mma_deallocate(TMP)
  end if
  ! Deleted orbitals:
  if (ND > 0) then
    call DCOPY_(ND,[Zero],0,OCC(IOCC+1),1)
    IOCC = IOCC+ND
    call DCOPY_(NB*ND,CMO(ICMO+1),1,CNAT(ICMO+1),1)
    ICMO = ICMO+NB*ND
  end if
end do

end subroutine NATORB_CASPT2
