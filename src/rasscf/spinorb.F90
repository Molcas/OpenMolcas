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

subroutine SPINORB(D,CMO,OCC)
! Purpose: diagonalize the spin density matrix (D) to
! obtain the eigenvectors (EVEC) and the eigenvalues (EVAL).
! Then the natural spinorbitals (CMONSO) are computed
! (only active).

use PrintLevel, only: DEBUG
use output_ras, only: LF, IPRLOC
use general_data, only: NSYM, NASH, NBAS, NFRO, NISH
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 D(*), CMO(*), OCC(*)
character(len=16) :: ROUTINE = 'SPINORB '
real*8, allocatable :: W1(:,:), W2(:,:)
integer :: I, IPCMO, IPDEN, IPOCC, iPrLev, iSym, NA, NB, NF, NI
integer :: IDIAG

! Local print level (if any)
IPRLEV = IPRLOC(6)
if (IPRLEV >= DEBUG) write(LF,*) ' Entering ',ROUTINE

IPDEN = 1
IPCMO = 1
IPOCC = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  NF = NFRO(ISYM)
  NI = NISH(ISYM)
  if (NB /= 0) then
    NA = NASH(ISYM)
    if (NA /= 0) then
      call mma_allocate(W1,NA,NA,Label='W1')
      call mma_allocate(W2,NB,NA,Label='W2')
      W1(:,:) = 0.0d0
      call DCOPY_(NA,[1.0d0],0,W1,NA+1)
      call Jacob(D(IPDEN),W1,NA,NA)
      IDIAG = 0
      do I=1,NA
        IDIAG = IDIAG+I
        OCC(IPOCC+NF+NI+I-1) = D(IPDEN+IDIAG-1)
      end do
      call DGEMM_('N','N',NB,NA,NA,1.0d0,CMO(IPCMO+(NF+NI)*NB),NB,W1,NA,0.0d0,W2,NB)
      call DCOPY_(NA*NB,W2,1,CMO(IPCMO+(NF+NI)*NB),1)
      call mma_deallocate(W2)
      call mma_deallocate(W1)
      IPDEN = IPDEN+NA*(NA+1)/2
    end if
    IPCMO = IPCMO+NB*NB
    IPOCC = IPOCC+NB
  end if
end do

end subroutine SPINORB
