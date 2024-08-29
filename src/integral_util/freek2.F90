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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine Freek2()
!***********************************************************************
!                                                                      *
!  Object: deallocate memory for pair entities                         *
!                                                                      *
! Called from: ClsSew                                                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. November '92                            *
!***********************************************************************

use k2_structure, only: free_k2data, Indk2, k2_Processed
use stdalloc, only: mma_deallocate

implicit none

if (allocated(Indk2)) then

  ! Deallocate k2 entities

  call Free_k2data()

  call mma_deallocate(Indk2)
  k2_processed = .false.

end if

return

end subroutine Freek2
