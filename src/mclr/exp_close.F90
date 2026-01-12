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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

subroutine Exp_Close()

use MCLR_Data, only: H0F, H0S, SBIDT
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(H0S,safe='*')
call mma_deallocate(H0F,safe='*')
call mma_deallocate(SBIDT,safe='*')

end subroutine Exp_Close
