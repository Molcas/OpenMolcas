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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nPLMO
real(kind=wp) :: PLMO(nPLMO)
integer(kind=iwp) :: mPLMO
logical(kind=iwp) :: Found
character(len=24) :: Label

Label = 'PLMO'
call qpg_dArray(Label,Found,mPLMO)
if ((.not. Found) .or. (mPLMO == 0)) call SysAbendMsg('get_plmo','Did not find:',Label)
if (nPLMO /= mPLMO) then
  write(u6,*) 'Get_PLMO: nPLMO/=mPLMO'
  write(u6,*) 'nPLMO=',nPLMO
  write(u6,*) 'mPLMO=',mPLMO
  call Abend()
end if
call get_dArray(Label,PLMO,nPLMO)

return

end subroutine Get_PLMO
