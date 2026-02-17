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

use Index_Functions, only: nTri_Elem
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NASH, NBAS, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: D(*), CMO(*), OCC(*)
integer(kind=iwp) :: I, IDIAG, IPCMO, IPDEN, IPOCC, iPrLev, iSym, NA, NB, NF, NI
real(kind=wp), allocatable :: W1(:,:), W2(:,:)

! Local print level (if any)
IPRLEV = IPRLOC(6)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering SPINORB'

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
      call unitmat(W1,NA)
      call Jacob(D(IPDEN),W1,NA,NA)
      IDIAG = 0
      do I=1,NA
        IDIAG = IDIAG+I
        OCC(IPOCC+NF+NI+I-1) = D(IPDEN+IDIAG-1)
      end do
      call DGEMM_('N','N',NB,NA,NA,One,CMO(IPCMO+(NF+NI)*NB),NB,W1,NA,Zero,W2,NB)
      CMO(IPCMO+(NF+NI)*NB:IPCMO+(NF+NI+NA)*NB-1) = pack(W2,.true.)
      call mma_deallocate(W2)
      call mma_deallocate(W1)
      IPDEN = IPDEN+nTri_Elem(NA)
    end if
    IPCMO = IPCMO+NB*NB
    IPOCC = IPOCC+NB
  end if
end do

end subroutine SPINORB
