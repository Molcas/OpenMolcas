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

subroutine f_inquire(filename,exists)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: filename
logical(kind=iwp), intent(out) :: exists
integer(kind=iwp) :: lRealName
character(len=256) :: RealName
#ifdef _SOLARIS_
integer(kind=iwp) :: n, irc
character(len=256) :: FTMP1, FTMP2

n = len(filename)
do
  if (filename(n:n) /= ' ') exit
  n = n-1
end do
n = n+1
ftmp1 = filename
ftmp1(n:n) = char(0)
exists = .true.
irc = -1
FTMP2 = ''
call PrgmTranslate(filename,RealName,lRealName)
FTMP1(1:lRealName) = RealName(1:lRealName)
!write(u6,*) 'before fndlnk', FTMP1
call fndlnk(irc,FTMP1,FTMP2)
if (irc /= 0) exists = .false.
#else

call PrgmTranslate(filename,RealName,lRealName)
!write(u6,*) 'Debug inquire:', RealName(1:lRealName)
inquire(file=RealName(1:lRealName),exist=exists)
#endif

return

end subroutine f_inquire
