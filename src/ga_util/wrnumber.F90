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

subroutine WrNumber(Name,Number)

use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: Number
character(len=*), intent(inout) :: Name
character(len=10) format
integer(kind=iwp) Limit, iTens
real(kind=wp) ANumber

format = ' '
if (Number >= 0) then
  Limit = 0
  do iTens=0,99
    Limit = Limit+9*(10**iTens)
    if (Number <= Limit) then
      write(format,'(A,I4.4,A)') '(I',iTens+1,')'
      write(Name,format) Number
      return
    end if
  end do
else
  ANumber = -dble(Number)
  Limit = 0
  do iTens=0,99
    Limit = Limit+9*(10**iTens)
    if (ANumber <= dble(Limit)) then
      write(format,'(A,I4.4,A)') '(A,I',iTens+1,')'
      write(Name,format) '-',ANumber
      return
    end if
  end do
end if
call SysAbendMsg('wrnumber','Number too large',' ')

end subroutine WrNumber
