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
      SubRoutine Cho_P_GetMQ(MQ,l_MQ,LstQSP,nQSP)

      Implicit None
      Integer l_MQ, nQSP
      Real*8  MQ(l_MQ)
      Integer LstQSP(nQSP)
#include "cho_para_info.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_P_GetMQ')

C --- In parallel:
C --- This code only works if MxShpr is set to 1
C ---      otherwise each node reads a slice
C ---      of MQ and thus a more sophisticated
C ---      "synchronization" would be needed
C ----------------------------------------------

      If (Cho_Real_Par) Then
         If (nQSP .gt. 1) Then
            Call Cho_Quit('Oops! Bug detected in '//SecNam,103)
         End If
         Call Cho_dZero(MQ,l_MQ)
         Call Cho_p_QualSwp()
         Call Cho_GetMQ(MQ,l_MQ,LstQSP,nQSP)
         Call Cho_p_QualSwp()
         Call Cho_GAdGop(MQ,l_MQ,'+')
      Else
         Call Cho_GetMQ(MQ,l_MQ,LstQSP,nQSP)
      End If

      End
