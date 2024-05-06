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

module wadr

! These arrays are used for the SXCTL part of the code:
!   BM, DIA, DIAF, F1, F2, SXG, SXH, SXN
! These arrays are used for the TRACTL2 part of the code:
!   CMO, D1A, D1I, FA, FI, OccN

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: NLX, nPWXY
real(kind=wp), allocatable :: BM(:), CMO(:), D1A(:), D1I(:), DIA(:), DIAF(:), DMAT(:), DSPN(:), F1(:), F2(:), FA(:), FI(:), &
                              FMO(:), FockOcc(:), OccN(:), PA(:), PMAT(:), SXG(:), SXH(:), SXN(:), TUVX(:)

public :: BM, CMO, D1A, D1I, DIA, DIAF, DMAT, DSPN, F1, F2, FA, FI, FMO, FockOcc, NLX, nPWXY, OccN, PA, PMAT, SXG, SXH, SXN, TUVX

end module wadr
