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

subroutine Get_DLAO(DLAO,nDLAO)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nDLAO
real(kind=wp) :: DLAO(nDLAO)
integer(kind=iwp) :: mDLAO
logical(kind=iwp) :: Found
character(len=24) :: Label

Label = 'DLAO'
call qpg_dArray(Label,Found,mDLAO)
if ((.not. Found) .or. (mDLAO == 0)) call SysAbendMsg('get_dlao','Did not find:',Label)
if (nDLAO /= mDLAO) then
  write(u6,*) 'Get_DLAO: nDLAO/=mDLAO'
  write(u6,*) 'nDLAO=',nDLAO
  write(u6,*) 'mDLAO=',mDLAO
  call Abend()
end if

call Get_dArray(Label,DLAO,nDLAO)

return

end subroutine Get_DLAO
