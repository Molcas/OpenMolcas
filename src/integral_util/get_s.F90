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

subroutine Get_S(icol,str,n)

use getline_mod, only: iend, istrt, line, ncol
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCol, n
character(len=*), intent(out) :: str(n)
integer(kind=iwp) :: i, ic

ic = icol
do i=1,n
  if (ic <= ncol) then
    if (iend(ic) >= istrt(ic)) then
      str(i) = line(istrt(ic):iend(ic))
    else
      str(i) = ' '
    end if
    ic = ic+1
  else
    write(u6,210) icol+n-1,line
    call FindErrorLine()
    call WarningMessage(2,'Error in Get_S')
    call Quit_OnUserError()
  end if
end do

210 format(/' ERROR IN GET_S: TRYING TO READ',i4,' STRINGS'/1x,a)

end subroutine Get_S
