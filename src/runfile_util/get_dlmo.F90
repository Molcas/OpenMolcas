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

subroutine Get_DLMO(DLMO,nDLMO)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nDLMO
real(kind=wp) :: DLMO(nDLMO)
integer(kind=iwp) :: mDLMO
logical(kind=iwp) :: Found
character(len=24) :: Label

Label = 'DLMO'
call qpg_dArray(Label,Found,mDLMO)
if ((.not. Found) .or. (mDLMO == 0)) call SysAbendMsg('get_dlmo','Did not find:',Label)
if (mDLMO /= nDLMO) then
  write(u6,*) 'Get_DLMO: nDLMO/=mDLMO'
  write(u6,*) 'nDLMO=',nDLMO
  write(u6,*) 'mDLMO=',mDLMO
  call Abend()
end if
call Get_dArray(Label,DLMO,nDLMO)

return

end subroutine Get_DLMO
