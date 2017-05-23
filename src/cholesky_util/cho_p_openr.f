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
      SubRoutine Cho_P_OpenR(iOpt)
C
C     Purpose: open (iOpt=1) or close (iOpt=2) files for storing global
C              reduced set indices.
C
      Implicit None
      Integer iOpt
#include "choglob.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_P_OpenR')

      Character*5 FNRed

      If (iOpt .eq. 1) Then
         LuRed_G = 7
         FNRed = 'CHRED'
         Call DAName_MF_WA(LuRed_G,FNRed)
      Else If (iOpt .eq. 2) Then
         If (LuRed_G .gt. 0) Then
            Call DAClos(LuRed_G)
         End If
      Else
         Call Cho_Quit('iOpt error in '//SecNam,104)
      End If

      End
