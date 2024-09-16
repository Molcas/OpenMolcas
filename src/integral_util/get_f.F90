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

subroutine Get_F(icol,val,n)

use getline_mod, only: iEnd, iStrt, Line, nCol
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCol, n
real(kind=wp), intent(out) :: val(n)
integer(kind=iwp) :: i, i1, i2, ic, istatus
character(len=80) :: string

ic = icol
do i=1,n
  if (ic <= ncol) then
    i1 = istrt(ic)
    i2 = iend(ic)
    if (i1 <= i2) then
      string = ' '
      string(len(string)+i1-i2:) = line(i1:i2)
      read(string,'(F80.0)',iostat=istatus) val(i)
      if (istatus /= 0) exit
    else
      val(i) = Zero
    end if
    ic = ic+1
  else
    write(u6,110) icol+n-1,line
    exit
  end if
end do

if (i <= n) then
  call FindErrorLine()
  call WarningMessage(2,'Error in Get_F')
  call Quit_OnUserError()
end if

return

110 format(/' ERROR IN GET_F: TRYING TO READ',i4,' VALUES'/1x,a)

end subroutine Get_F
