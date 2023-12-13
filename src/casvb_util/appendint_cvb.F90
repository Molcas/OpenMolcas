!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine appendint_cvb(c,nmbr,iskip)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(inout) :: c
integer(kind=iwp), intent(in) :: nmbr, iskip
integer(kind=iwp) :: ibegin, iend, itens, limit, mnumber
character(len=10) :: frmt

ibegin = len_trim(c)+1+iskip
iend = len(c)

frmt = ' '
if (nmbr >= 0) then
  limit = 0
  do itens=0,99
    limit = limit+9*(10**itens)
    if (nmbr <= limit) then
      write(frmt,'(a,i4.4,a)') '(i',itens+1,')'
      write(c(ibegin:iend),frmt) nmbr
      return
    end if
  end do
else
  mnumber = -nmbr
  limit = 0
  do itens=0,99
    limit = limit+9*(10**itens)
    if (mnumber <= limit) then
      write(frmt,'(a,i4.4,a)') '(a,i',itens+1,')'
      write(c(ibegin:iend),frmt) '-',mnumber
      return
    end if
  end do
end if
write(u6,*) ' Number too large in appendint :',nmbr
call abend_cvb()

end subroutine appendint_cvb
