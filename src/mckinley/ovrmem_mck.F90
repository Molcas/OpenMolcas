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

subroutine OvrMem_mck( &
#                     define _CALLING_
#                     include "mem_interface.fh"
                     )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: lr_

#include "macros.fh"
unused_var(lr)

lr_ = 0
nHer = (la+lb+lr_+3)/2
Mem = 3*nHer*(la+2)+3*nHer*(lb+2)+3*nHer*(lr_+2)+3*(la+2)*(lb+2)*(lr_+1)+2+nTri_Elem1(la)*nTri_Elem1(lb)*2

return

end subroutine OvrMem_mck
