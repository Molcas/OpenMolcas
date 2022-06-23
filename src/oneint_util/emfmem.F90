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

subroutine EMFMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

nHer = (la+lb+lr+2)/2
Mem = 3*nHer*(la+1+lr)*2+3*nHer*(lb+1+lr)*2+3*(la+1+lr)*(lb+1+lr)*2
if (lr == 1) then
  Mem = Mem+6*(la+1)*(lb+1)*2+2+nElem(la)*nElem(lb)*nElem(lr)*12
else
  Mem = Mem+nElem(la)*nElem(lb)*nElem(lr)*2
end if

return

end subroutine EMFMem
