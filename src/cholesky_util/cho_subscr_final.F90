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

subroutine Cho_SubScr_Final()
!
! Purpose: finalize (de-allocate memory) screening in vector
!          subtraction.

use Cholesky, only: DSPNm, DSubScr
use stdalloc, only: mma_deallocate

implicit none

if (allocated(DSPNm)) call mma_deallocate(DSPNm)

if (allocated(DSubScr)) call mma_deallocate(DSubScr)

end subroutine Cho_SubScr_Final
