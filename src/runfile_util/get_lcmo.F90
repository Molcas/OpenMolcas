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

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nLCMO
real(kind=wp) :: LCMO(nLCMO)
integer(kind=iwp) :: mLCMO
logical(kind=iwp) :: Found
character(len=24) :: Label

Label = 'LCMO'
call qpg_dArray(Label,Found,mLCMO)
if ((.not. Found) .or. (mLCMO == 0)) call SysAbendMsg('get_lcmo','Did not find:',Label)
if (nLCMO /= mLCMO) then
  write(u6,*) 'Get_LCMO: nLCMO/=mLCMO'
  write(u6,*) 'nLCMO=',nLCMO
  write(u6,*) 'mLCMO=',mLCMO
  call Abend()
end if
call get_dArray(Label,LCMO,nLCMO)

return

end subroutine Get_LCMO
