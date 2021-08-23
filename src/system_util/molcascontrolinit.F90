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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************

subroutine MolcasControlInit(label)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Label
character(len=512) :: tmp
character(len=32) :: My
integer(kind=iwp) :: i, ic, islast, Lu
integer(kind=iwp), parameter :: nLines = 20
character(len=*), parameter :: filename = 'molcas.control'
integer(kind=iwp), external :: StrnLn

iC = 0
tmp = Label(1:len(Label))
Lu = 1
call molcas_open(Lu,filename)
write(Lu,'(a)') '# Molcas control file: change # to ! to activate.'
islast = 0

do
  i = index(tmp,',')
  My = ' '
  if (i > 0) then
    My(1:i-1) = tmp(1:i-1)
    tmp = tmp(i+1:)
  else
    My = trim(tmp)
    islast = 1
  end if
  iC = iC+1
  if (ic > nLines) call abend()
  i = StrnLn(My)
  if (index(My,'=') == 0) then
    i = i+1
    My(i:i) = '='
  end if
  write(Lu,'(a,a)') '#',My(1:i)
  if (islast /= 0) exit
end do
close(Lu)

return

end subroutine MolcasControlInit
