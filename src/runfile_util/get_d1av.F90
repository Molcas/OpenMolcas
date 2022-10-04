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

subroutine Get_D1AV(D1AV,nD1AV)

implicit real*8(A-H,O-Z)
character*24 Label
logical Found
real*8 D1AV(nD1AV)

Label = 'D1av'
call qpg_dArray(Label,Found,mD1AV)
if ((.not. Found) .or. (mD1AV == 0)) call SysAbendMsg('Get_D1AV','Did not find:',Label)
if (nD1AV /= mD1AV) then
  write(6,*) 'Get_D1AV: nD1AV/=mD1AV'
  write(6,*) 'nD1AV=',nD1AV
  write(6,*) 'mD1AV=',mD1AV
  call Abend()
end if

call Get_dArray(Label,D1AV,nD1AV)

return

end subroutine Get_D1AV
