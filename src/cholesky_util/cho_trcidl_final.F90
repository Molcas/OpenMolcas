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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_TrcIdl_Final()
!
! Thomas Bondo Pedersen, May 2010.
!
! Deallocate array for tracing idle processors

use Cholesky, only: Idle
use stdalloc, only: mma_deallocate

implicit none

if (allocated(Idle)) call mma_deallocate(Idle)

end subroutine Cho_TrcIdl_Final
