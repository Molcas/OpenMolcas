************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine DecideOnDF(DoDF)
      Implicit None
      Logical DoDF

      Integer iOption

      Call Get_iScalar('System BitSwitch',iOption)
      DoDF=IAND(iOption,1024).eq.1024

#if defined (_DEBUG_)
      Write(6,*) '>>> Exit from DecideOnDF:'
      Write(6,*) '    System Bit Switch = ',iOption
      Write(6,*) '    DoDF = ',DoDF
#endif

      End
