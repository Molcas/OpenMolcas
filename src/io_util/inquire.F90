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

subroutine f_inquire(NAME,exist)
character*(*) name
logical exist
character*256 RealName
#ifdef _SOLARIS_
character*256 FTMP1, FTMP2
integer n, irc

n = len(name)
200 continue
if (name(n:n) == ' ') then
  n = n-1
  Go To 200
end if
n = n+1
ftmp1 = name
ftmp1(n:n) = char(0)
exist = .true.
irc = -1
FTMP2 = ''
call PrgmTranslate(Name,RealName,lRealName)
FTMP1(1:lRealName) = RealNAME(1:lRealName)
!write(6,*) 'before fndlnk', FTMP1
call fndlnk(irc,FTMP1,FTMP2)
if (irc /= 0) exist = .false.
#else

call PrgmTranslate(Name,RealName,lRealName)
!write(6,*) 'Debug inquire:', RealNAME(1:lRealName)
inquire(file=RealNAME(1:lRealName),exist=exist)
#endif

return

end subroutine f_inquire
