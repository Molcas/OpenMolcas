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

subroutine OMQMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

#include "mem_interface.fh"
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

call MltMmP(nOrder,MemP,la,lb+1,lr-1)
! Not sure what this does. nHer is set on the outside anyway???
nHer = nOrder
! not same order (1 not 3) hence lr-2
call MltMmP(nOrder,MemN,la,lb,lr-2)
Mem = max(MemP,MemN)
if (lb > 0) then
  call MltMmP(nOrder,MemM,la,lb-1,lr-1)
  ! For L - 1 Component
  Mem = max(Mem,MemM)+nElem(la)*nElem(lb-1)*6
end if
! For dipole term ( L + 0 Component)
Mem = Mem+nElem(la)*nElem(lb)*3
! For L + 1 Component
Mem = Mem+1+nElem(la)*nElem(lb+1)*6
Mem = Mem+nElem(la)*nElem(lb)*9 ! final term

return

end subroutine OMQMem
