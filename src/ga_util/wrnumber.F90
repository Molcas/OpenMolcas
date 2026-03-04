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

subroutine WrNumber(Nam,Nmbr)

use Definitions, only: wp, iwp

implicit none
character(len=*), intent(inout) :: Nam
integer(kind=iwp), intent(in) :: Nmbr
integer(kind=iwp) :: iTens, Limit
real(kind=wp) :: ANumber
character(len=10) :: frmt

frmt = ' '
if (Nmbr >= 0) then
  Limit = 0
  do iTens=0,99
    Limit = Limit+9*(10**iTens)
    if (Nmbr <= Limit) then
      write(frmt,'(A,I4.4,A)') '(I',iTens+1,')'
      write(Nam,frmt) Nmbr
      return
    end if
  end do
else
  ANumber = real(-Nmbr,kind=wp)
  Limit = 0
  do iTens=0,99
    Limit = Limit+9*(10**iTens)
    if (ANumber <= real(Limit,kind=wp)) then
      write(frmt,'(A,I4.4,A)') '(A,I',iTens+1,')'
      write(Nam,frmt) '-',ANumber
      return
    end if
  end do
end if
call SysAbendMsg('wrnumber','Number too large',' ')

end subroutine WrNumber
