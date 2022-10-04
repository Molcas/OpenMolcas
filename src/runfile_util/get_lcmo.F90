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

subroutine Get_LCMO(LCMO,nLCMO)

implicit real*8(A-H,O-Z)
character*24 Label
logical Found
real*8 LCMO(nLCMO)

Label = 'LCMO'
call qpg_dArray(Label,Found,mLCMO)
if ((.not. Found) .or. (mLCMO == 0)) call SysAbendMsg('get_lcmo','Did not find:',Label)
if (nLCMO /= mLCMO) then
  write(6,*) 'Get_LCMO: nLCMO/=mLCMO'
  write(6,*) 'nLCMO=',nLCMO
  write(6,*) 'mLCMO=',mLCMO
  call Abend()
end if
call get_dArray(Label,LCMO,nLCMO)

return

end subroutine Get_LCMO
