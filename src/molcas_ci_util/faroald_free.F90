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

subroutine faroald_free()
! The finalization subroutine lives outside of the
! faroald module so that it can be called separately.

use faroald, only: ex1_a, ex1_b, mma_deallocate

implicit none

if (allocated(ex1_a)) call mma_deallocate(ex1_a)
if (allocated(ex1_b)) call mma_deallocate(ex1_b)

end subroutine faroald_free
