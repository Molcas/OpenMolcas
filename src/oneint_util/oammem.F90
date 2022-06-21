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

subroutine OAMMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

#include "mem_interface.fh"
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

call MltMmP(nOrder,Mem,la,lb+1,lr-1)
nHer = nOrder
if (lb > 0) then
  call MltMmP(nOrder,MmMltP,la,lb-1,lr-1)
  Mem = max(Mem,MmMltP)+nElem(la)*nElem(lb-1)*3
end if
Mem = Mem+1+nElem(la)*nElem(lb+1)*3
Mem = Mem+nElem(la)*nElem(lb)*3

return

end subroutine OAMMem
