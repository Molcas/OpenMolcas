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
!*****************************************************************************
!                                                                            *
! Author:   Valera Veryazov 2000-2016                                        *
!           Theoretical Chemistry                                            *
!           Lund University                                                  *
!           Sweden                                                           *
!                                                                            *
!*****************************************************************************

subroutine StdIn_Name(Name)

character*(*) Name
character*132 Line

nName = len(Name)
if (nName /= 16) then
  write(6,*) 'StdIn_Name: Wrong length of character Name'
  call Abend()
end if
!write(6,*) 'nName=',nName

Name = 'Stdin.  '
call GetEnvf('EMIL_RC2',Line)
read(Line,'(I132.132)') Iter
!write(6,*) 'Line=',Line
!write(6,*) 'Iter=',Iter
Iter = Iter+1
if (Line(1:1) == ' ') then
  Name(7:7) = '2'
else if (Iter <= 9) then
  write(Name(7:7),'(I1)') Iter
else if (Iter <= 99) then
  write(Name(7:8),'(I2)') Iter
else
  write(6,*) 'StdIn_Name: Error in Line!'
  call Abend()
end if
Line = ' '
call GetEnvf('EMIL_InLoop',Line)
ib = -1
ie = -1
i = 1
10 continue
if (Line(i:i) /= ' ' .and. ib == -1) then
  ib = i
  goto 20
end if
if (Line(i:i) == ' ' .and. ib > 0) then
  ie = i
  goto 30
end if
20 i = i+1
goto 10
30 Name(index(Name,' '):) = '.'//Line(ib:ie)
!write(6,*) '>',Line,'<',ib,ie
!write(6,*) 'StdIn=',Name

return

end subroutine StdIn_Name
