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
      Subroutine LDF_AddConstraint(Constraint)
      Implicit None
      Integer Constraint
#include "localdf.fh"
      If (Constraint.le.-1) Then
         LDF_Constraint=-1
      Else If (Constraint.le.LDF_MXCONSTRAINT) Then
         LDF_Constraint=Constraint
      Else
         Call WarningMessage(2,'LDF constraint not recognized')
         Write(6,'(A,I10,A,I10)')
     &   'Constraint=',Constraint,' > ',LDF_MXCONSTRAINT
         Call Quit_OnUserError()
      End If
      End
