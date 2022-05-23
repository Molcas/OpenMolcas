************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000-2016, Valera Veryazov                             *
************************************************************************
      subroutine f_inquire(NAME,exist)
      Character*(*) name
      Logical exist
      Character*256 RealName
#ifdef _SOLARIS_
      Character*256 FTMP1,FTMP2
      Integer n,irc
      n=Len(name)
200   If(name(n:n).eq.' ') Then
         n=n-1
         Go To 200
      End If
      n=n+1
      ftmp1=name
      ftmp1(n:n)=Char(0)
      exist=.true.
      irc=-1
      FTMP2='                           '//
     &      '                           '
      Call PrgmTranslate(Name,RealName,lRealName)
      FTMP1(1:lRealName)=RealNAME(1:lRealName)
c         print *,'before fndlnk', FTMP1
      call fndlnk(irc,FTMP1,FTMP2)
      if (irc.ne.0) exist=.false.
#else

        Call PrgmTranslate(Name,RealName,lRealName)
c      print *,'Debug inquire:', RealNAME(1:lRealName)
        inquire(file=RealNAME(1:lRealName),exist=exist)
#endif
      return
      end
