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

subroutine Get_PLMO(PLMO,nPLMO)

implicit real*8(A-H,O-Z)
character*24 Label
logical Found
real*8 PLMO(nPLMO)

Label = 'PLMO'
call qpg_dArray(Label,Found,mPLMO)
if ((.not. Found) .or. (mPLMO == 0)) call SysAbendMsg('get_plmo','Did not find:',Label)
if (nPLMO /= mPLMO) then
  write(6,*) 'Get_PLMO: nPLMO/=mPLMO'
  write(6,*) 'nPLMO=',nPLMO
  write(6,*) 'mPLMO=',mPLMO
  call Abend()
end if
call get_dArray(Label,PLMO,nPLMO)

return

end subroutine Get_PLMO
