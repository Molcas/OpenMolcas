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

subroutine Get_P2MOt(P2MO,nP2MO)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
character*24 Label
logical Found
real*8 P2MO(nP2MO)

Label = 'P2MOT'
call qpg_dArray(Label,Found,mP2MO)
if ((.not. Found) .or. (mP2MO == 0)) call SysAbendmsg('Get_P2MOt','Did not find:',label)
if (nP2MO /= mP2MO) then
  write(6,*) 'Get_P2MO: nP2MO/=mP2MO'
  write(6,*) 'mP2MO=',mP2MO
  write(6,*) 'nP2MO=',nP2MO
  call Abend()
end if
call get_dArray(Label,P2MO,nP2MO)

return

end subroutine Get_P2MOt
