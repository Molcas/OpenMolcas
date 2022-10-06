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

subroutine Get_dExcdRa_X(dExcdRa,ndExcdRa)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ndExcdRa
real(kind=wp), intent(out) :: dExcdRa(ndExcdRa)
integer(kind=iwp) :: mdExcdRa
logical(kind=iwp) :: Found
character(len=24) :: Label

Label = 'dExcdRa'
call qpg_dArray(Label,Found,mdExcdRa)
if ((.not. Found) .or. (mdExcdRa == 0)) call SysAbendmsg('Get_dExcdRa_X','Did not find:',Label)
if (mdExcdRa /= ndExcdRa) then
  write(u6,*) 'mdExcdRa/=ndExcdRa'
  write(u6,*) mdExcdRa,'/=',ndExcdRa
  call AbEnd()
end if
call Get_dArray(Label,dExcdRa,ndExcdRa)

return

end subroutine Get_dExcdRa_X
