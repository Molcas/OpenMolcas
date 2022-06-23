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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine MVeMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"

nHer = (la+lb+lr+2)/2
Mem = 3*nHer*(la+3)+3*nHer*(lb+3)+3*nHer*(lr-3)+3*(la+3)*(lb+3)*(lr-3)+3*(la+1)*(lb+1)*2+3*(la+1)*(lb+1)+1+1

return

end subroutine MVeMem
