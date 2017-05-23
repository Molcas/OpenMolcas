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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Quit(SecNam,Str1,Str2)
C
C     Thomas Bondo Pedersen, Jan. 2005.
C
C     Purpose: stop execution using ChoMP2_Quit and print the trace
C              stack.
C
      Implicit None
      Character*(*) SecNam, Str1, Str2

      Call qTrace
      Call SysAbendMsg(SecNam,Str1,Str2)

      End
