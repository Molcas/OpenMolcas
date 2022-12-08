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

subroutine RFMmG( &
#                define _CALLING_
#                include "mem_interface.fh"
                )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"

nHer = (la+lb+lr+1+2)/2
Mem = 3*nHer*(la+2)+ &
      3*nHer*(lb+2)+ &
      3*nHer*(lr+1)+ &
      3*(la+2)*(lb+2)*(lr+1)+4+3*nHer

return

end subroutine RFMmG
