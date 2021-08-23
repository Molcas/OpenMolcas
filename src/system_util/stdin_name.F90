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
!***********************************************************************
!                                                                      *
! Author:   Valera Veryazov 2000-2016                                  *
!           Theoretical Chemistry                                      *
!           Lund University                                            *
!           Sweden                                                     *
!                                                                      *
!***********************************************************************

subroutine StdIn_Name(Name)

use Definitions, only: iwp, u6

implicit none
character(len=16), intent(out) :: Name
character(len=132) :: Line
integer(kind=iwp) :: i, ib, ie, Iter, nName

nName = len(Name)
if (nName /= 16) then
  write(u6,*) 'StdIn_Name: Wrong length of character Name'
  call Abend()
end if
!write(u6,*) 'nName=',nName

Name = 'Stdin.  '
call GetEnvf('EMIL_RC2',Line)
read(Line,'(I132.132)') Iter
!write(u6,*) 'Line=',Line
!write(u6,*) 'Iter=',Iter
Iter = Iter+1
if (Line(1:1) == ' ') then
  Name(7:7) = '2'
else if (Iter <= 9) then
  write(Name(7:7),'(I1)') Iter
else if (Iter <= 99) then
  write(Name(7:8),'(I2)') Iter
else
  write(u6,*) 'StdIn_Name: Error in Line!'
  call Abend()
end if
Line = ' '
call GetEnvf('EMIL_InLoop',Line)
ib = -1
ie = -1
do i=1,len(Line)
  if ((Line(i:i) /= ' ') .and. (ib == -1)) then
    ib = i
  else if ((Line(i:i) == ' ') .and. (ib > 0)) then
    ie = i
    exit
  end if
end do
Name(index(Name,' '):) = '.'//Line(ib:ie)
!write(u6,*) '>',Line,'<',ib,ie
!write(u6,*) 'StdIn=',Name

return

end subroutine StdIn_Name
