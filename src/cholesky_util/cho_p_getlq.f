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
      SubRoutine Cho_P_GetLQ(QVec,l_QVec,LstQSP,nQSP)

      Implicit None
      Integer l_QVec, nQSP
      Real*8, Target::  QVec(l_Qvec)
      Integer LstQSP(nQSP)
#include "cho_para_info.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_P_GetLQ')

C --- In parallel:
C --- This code only works if MxShpr is set to 1
C ---      otherwise each node computes a slice
C ---      of QVec and thus a more sophisticated
C ---      "synchronization" would be needed

      If (Cho_Real_Par) Then
         If (nQSP .gt. 1) Then
            Call Cho_Quit('Oops! Bug detected in '//SecNam,103)
         End If
         Call Cho_dZero(QVec,l_Qvec)
         Call Cho_p_QualSwp()
         Call Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
         Call Cho_p_QualSwp()
         Call Cho_GAdGOp(QVec,l_QVec,'+') ! sync. array
      Else
         Call Cho_GetLQ(QVec,l_QVec,LstQSP,nQSP)
      End If

      End
