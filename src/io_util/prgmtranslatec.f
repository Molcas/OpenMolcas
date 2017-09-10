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

* Wrapper for the OpenMolcas PrgmTranslate_Mod routine, so it works
* without using the module (i.e., from C). Some C-to-Fortran conversion
* needs to be done, which feels quite hackish.

#ifndef _HAVE_EXTRA_
#define MAXSTR 1024
      Subroutine PrgmTranslateC(InStr,l1,OutStr,l2,Par)
      Use ISO_C_Binding, Only: C_CHAR, C_NULL_CHAR
      Use Prgm
      Implicit None
      Character(Kind=C_CHAR), Intent(In) :: InStr(*)
      Character(Kind=C_CHAR), Intent(Out) :: OutStr(*)
      Integer, Intent(In) :: Par, l1
      Integer, Intent(Out) :: l2
      Character (Len=MAXSTR) :: TmpStr, TmpStr2
      Integer :: i
      TmpStr=''
      Do i=1,l1
        TmpStr(i:i)=InStr(i)
      End Do
      Call PrgmTranslate_Mod(TmpStr,l1,TmpStr2,l2,Par)
      Do i=1,l2
        OutStr(i)=TmpStr2(i:i)
      End Do
      OutStr(l2+1:l2+1)=C_NULL_CHAR
      End Subroutine PrgmTranslateC
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_prgmtranslatec
      End
#endif
