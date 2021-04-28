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

subroutine PPMem( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Index_util, only: nTri0Elem
use Definitions, only: iwp

implicit none
#define _USE_WP_
#include "mem_interface.fh"
integer(kind=iwp) :: intmax

nHer = 0
Mem = 0
intmax = max(nTri0Elem(la),nTri0Elem(lb))
intmax = intmax**2
Mem = Mem+3*intmax

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lr)

end subroutine PPMem
