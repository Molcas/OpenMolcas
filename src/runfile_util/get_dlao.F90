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

implicit real*8(A-H,O-Z)
character*24 Label
logical Found
real*8 DLAO(nDLAO)

Label = 'DLAO'
call qpg_dArray(Label,Found,mDLAO)
if ((.not. Found) .or. (mDLAO == 0)) call SysAbendMsg('get_dlao','Did not find:',Label)
if (nDLAO /= mDLAO) then
  write(6,*) 'Get_DLAO: nDLAO/=mDLAO'
  write(6,*) 'nDLAO=',nDLAO
  write(6,*) 'mDLAO=',mDLAO
  call Abend()
end if

call Get_dArray(Label,DLAO,nDLAO)

return

end subroutine Get_DLAO
