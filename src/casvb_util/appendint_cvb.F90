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

subroutine appendint_cvb(c,number,iskip)

implicit real*8(a-h,o-z)
character*(*) c
character*10 format

ibegin = len_trim_cvb(c)+1+iskip
iend = len(c)

format = ' '
if (number >= 0) then
  limit = 0
  do itens=0,99
    limit = limit+9*(10**itens)
    if (number <= limit) then
      write(format,'(a,i4.4,a)') '(i',itens+1,')'
      write(c(ibegin:iend),format) number
      return
    end if
  end do
else
  mnumber = -number
  limit = 0
  do itens=0,99
    limit = limit+9*(10**itens)
    if (mnumber <= limit) then
      write(format,'(a,i4.4,a)') '(a,i',itens+1,')'
      write(c(ibegin:iend),format) '-',mnumber
      return
    end if
  end do
end if
write(6,*) ' Number too large in appendint :',number
call abend_cvb()

end subroutine appendint_cvb
