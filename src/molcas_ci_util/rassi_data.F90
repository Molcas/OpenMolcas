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

module rassi_data

! NISHT  - NR OF  INACTIVE ORBITALS, TOTAL
! NASHT  - NR OF    ACTIVE ORBITALS, TOTAL
! NOSHT  - NR OF  OCCUPIED ORBITALS, TOTAL
! NSSHT  - NR OF SECONDARY ORBITALS, TOTAL
! NISH(8),...,NSSH(8), AS ABOVE, BUT BY SYMMETRY TYPE.
! NBST, NBAS(8) - SIMILAR, NR OF BASIS FUNCTIONS.
! NBMX   - MAX NR OF BASIS FUNCTIONS OF ANY SPECIFIC SYMMETRY.
! NBTRI  - TOTAL SIZE OF TRIANGULAR SYMMETRY BLOCKS OF BASIS FNCS.
! NBSQ   - D:O, SQUARE SYMMETRY BLOCKS.
! NBSQPR - ACCUMULATED NR OF SQUARE SYMMETRY BLOCKS OF PREVIOUS
!          SYMMETRY TYPES. NBSQPR(1)=0.
! NCMO   - SIZE OF CMO COEFFICIENT ARRAYS, = SUM(NOSH(I)*NBAS(I)).
! NTDMAB - SIZE OF TRANS.D. MATRIX IN BIORTHONORMAL MO BASIS.
! NTDMZZ - SIZE OF TRANS.D. MATRIX IN AO BASIS.
! NSXY   - SIZE OF MO OVERLAP ARRAY.
! NTRA   - SIZE OF TRANSFORMATION COEFFICIENT ARRAY.
! ChFracMem - fraction of memory for the Cholesky vectors buffer

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: NAES(8), NASH(8), NASHT, NBASF(8), NBMX, NBSQ, NBSQPR(8), NBST, NBTRI, NCMO, NDEL(8), NFRO(8), NISH(8), &
                     NISHT, NOSH(8), NOSHT, NSSH(8), NSSHT, NSXY, NTDMAB, NTDMZZ, NTRA
character(len=8) :: WFTYPE
real(kind=wp) :: ChFracMem, ENUC

public :: ChFracMem, ENUC, NAES, NASH, NASHT, NBASF, NBMX, NBSQ, NBSQPR, NBST, NBTRI, NCMO, NDEL, NFRO, NISH, NISHT, NOSH, NOSHT, &
          NSSH, NSSHT, NSXY, NTDMAB, NTDMZZ, NTRA, WFTYPE

end module rassi_data
