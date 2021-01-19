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
      SubRoutine Cho_P_SetAddr()
C
C     Purpose: set initial disk adresses for local as well as global
C              reduced sets.
C
      use ChoSwp, only: InfRed, InfRed_G
      use ChoSwp, only: InfVec, InfVec_G
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "cho_para_info.fh"
#include "choglob.fh"

      Character*13 SecNam
      Parameter (SecNam = 'Cho_P_SetAddr')

      Integer irc

      If (Cho_Real_Par) Then

C        The variable XnPass must be zero: restart is not possible.
C        ----------------------------------------------------------

         If (XnPass .ne. 0) Then
            Call Cho_Quit('XnPass>0 error in '//SecNam,104)
         End If

C        Global.
C        -------

         Call Cho_P_SetAddr_2(InfRed_G,InfVec_G,
     &                        MaxRed,MaxVec,SIZE(InfVec,2),nSym,irc)
         If (irc .ne. 0) Then
            Write(Lupri,*) SecNam,': Cho_P_SetAddr_2 returned ',irc
            Call Cho_Quit('Error in '//SecNam,104)
         End If

      End If

C     Local.
C     ------

      Call Cho_SetAddr(InfRed,InfVec,
     &                 MaxRed,MaxVec,SIZE(InfVec,2),nSym)

      End
************************************************************************
************************************************************************
************************************************************************
      SubRoutine Cho_P_SetAddr_2(InfRed,InfVec,MaxRed,MaxVec,N2,nSym,
     &                           irc)
      Implicit None
      Integer MaxRed, MaxVec, N2, nSym, irc
      Integer InfRed(MaxRed), InfVec(MaxVec,N2,nSym)

      Integer iSym

      irc = 0

      If (MaxRed .gt. 0) Then
         InfRed(1) = 0
      Else
         irc = 1
         Return
      End If

      If (MaxVec.gt.0 .and. N2.ge.4) Then
         Do iSym = 1,nSym
            InfVec(1,3,iSym) = 0
            InfVec(1,4,iSym) = 0
         End Do
      Else
         irc = 2
         Return
      End If

      End
