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

subroutine Free_DeDe_Funi()

use k2_arrays, only: DeDe, ipOffD
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(ipOffD)
call mma_deallocate(DeDe)

return

end subroutine Free_DeDe_Funi
