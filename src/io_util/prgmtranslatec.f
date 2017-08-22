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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************

* Wrapper for the OpenMolcas PrgmTranslateC routine, so it works without using the module

#ifndef _HAVE_EXTRA_
      Subroutine PrgmTranslateC(InStr,l1,OutStr,l2,Par)
      Use Prgm
      Implicit None
      Character (Len=*), Intent(In) :: InStr
      Character (Len=*), Intent(Out) :: OutStr
      Character (Len=Len(InStr)) :: TmpStr
      Integer, Intent(In) :: Par, l1
      Integer, Intent(Out) :: l2
      TmpStr=InStr
      If (Len(InStr) .gt. l1) TmpStr(l1+1:)=''
      Call PrgmTranslate_Mod(TmpStr,l1,OutStr,l2,Par)
      If (Len(OutStr) .gt. l2) OutStr(l2+1:l2+1)=Char(0)
      End Subroutine PrgmTranslateC
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_prgmtranslatec
      End
#endif
