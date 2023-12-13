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

subroutine KnEMem_GIAO( &
#                      define _CALLING_
#                      include "mem_interface.fh"
                      )

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
#include "mem_interface.fh"

nHer = ((la+1)+(lb+1)+lr+3)/2

Mem = 3*nHer*(la+2)+3*nHer*(lb+2)+3*nHer*(lr+2)+3*(la+2)*(lb+2)*(lr+2)+3*(la+1)*(lb+1)*(lr+2)+3*(la+1)*(lb+1)*2+1+1+ &
      nTri_Elem1(la)*nTri_Elem1(lb)*3

return

end subroutine KnEMem_GIAO
