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
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************

subroutine TRACHOSZ_FREE()

use CHOVEC_IO, only: NVLOC_CHOBATCH, IDLOC_CHOGROUP, NVGLB_CHOBATCH, IDGLB_CHOGROUP
use stdalloc, only: mma_deallocate

implicit none

call MMA_DEALLOCATE(NVLOC_CHOBATCH)
call MMA_DEALLOCATE(IDLOC_CHOGROUP)
call MMA_DEALLOCATE(NVGLB_CHOBATCH)
call MMA_DEALLOCATE(IDGLB_CHOGROUP)

end subroutine TRACHOSZ_FREE
