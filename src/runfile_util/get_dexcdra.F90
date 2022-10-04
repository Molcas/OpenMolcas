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

subroutine Get_dExcdRa(ipdExcdRa,ndExcdRa)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
character(LEN=24) Label
logical Found

Label = 'dExcdRa'
call qpg_dArray(Label,Found,ndExcdRa)
if ((.not. Found) .or. (ndExcdRa == 0)) call SysAbendmsg('Get_dExcdRa','Did not find:',Label)
call GetMem('dExcdRa','Allo','Real',ipdExcdRa,ndExcdRa)
call Get_dArray(Label,Work(ipdExcdRa),ndExcdRa)

return

end subroutine Get_dExcdRa
