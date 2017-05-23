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
      Real*8 Function OrbPhase(A,nA)
      Implicit None
      Integer nA, i
      Real*8 A(nA)
*
      OrbPhase=0.0D0
      Do i = 1, nA
         OrbPhase = OrbPhase + A(i)*DBLE(i)
      End Do
      If (OrbPhase.lt.0.0D0) Then
         Call DSCal_(nA,-1.0D0,A,1)
         OrbPhase = -OrbPhase
      End If
*
      Return
      End
