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

subroutine Get_D2AV(D2AV,nD2AV)

implicit real*8(A-H,O-Z)
character*24 Label
logical Found
real*8 D2AV(nD2AV)

Label = 'D2av'
call qpg_dArray(Label,Found,mD2AV)
if ((.not. Found) .or. (mD2AV == 0)) call SysAbendMsg('get_d2av','Did not find',Label)
if (nD2AV /= mD2AV) then
  write(6,*) 'Get_D2AV: nD2AV/=mD2AV'
  write(6,*) 'nD2AV=',nD2AV
  write(6,*) 'mD2AV=',mD2AV
  call Abend()
end if
call get_dArray(Label,D2AV,nD2AV)

return

end subroutine Get_D2AV
