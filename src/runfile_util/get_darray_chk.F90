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

subroutine Get_dArray_chk(Label,rData,nData)

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: nData
real(kind=wp), intent(out) :: rData(nData)
integer(kind=iwp) :: mData
logical(kind=iwp) :: Found

call qpg_dArray(Label,Found,mData)
if ((.not. Found) .or. (mData == 0)) then
  call SysAbendMsg('Get_dArray_chk','Did not find:',Label)
else if (nData /= mData) then
  write(u6,*) 'Get_dArray_chk: nData /= mData'
  write(u6,*) 'nData=',nData
  write(u6,*) 'mData=',mData
  call Abend()
end if
call Get_dArray(Label,rData,nData)

return

end subroutine Get_dArray_chk
