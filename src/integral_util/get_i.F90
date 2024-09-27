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

subroutine Get_I(icol,ival,n)

use getline_mod, only: iEnd, iStrt, Line, nCol
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCol, n
integer(kind=iwp), intent(out) :: ival(n)
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
      read(string,'(i80)',iostat=istatus) ival(i)
      if (istatus /= 0) exit
    else
      ival(i) = 0
    end if
    ic = ic+1
  else
    write(u6,210) icol+n-1,line
    exit
  end if
end do

if (i <= n) then
  call FindErrorLine()
  call WarningMessage(2,'Error in Get_I')
  call Quit_OnUserError()
end if

return

210 format(/' ERROR IN GET_I: TRYING TO READ',i4,' VALUES'/1x,a)

end subroutine Get_I
