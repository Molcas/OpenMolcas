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

subroutine D1Mem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Sizes_of_Seward, only: S

implicit real*8(A-H,O-Z)
#include "mem_interface.fh"

nHer = S%mCentr
Mem = 3*(la+1)*nHer+3*(lb+1)*nHer

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lr)

end subroutine D1Mem
