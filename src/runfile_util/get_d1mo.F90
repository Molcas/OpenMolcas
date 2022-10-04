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

subroutine Get_D1MO(D1MO,nDens)

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
character*24 Label
logical Found
integer nDens
real*8 D1MO(nDens)

Label = 'D1mo'
call qpg_dArray(Label,Found,mDens)
if ((.not. Found) .or. (nDens == 0)) call SysAbendMsg('get_d1mo','Did not find:',Label)
if (mDens /= nDens) then
  write(6,*) 'Get_D1MO: mDens/=nDens'
  write(6,*) 'mDens=',mDens
  write(6,*) 'nDens=',nDens
  call Abend()
end if
call Get_dArray(Label,D1MO,nDens)

return

end subroutine Get_D1MO
