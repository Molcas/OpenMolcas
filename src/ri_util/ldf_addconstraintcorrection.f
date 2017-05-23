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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_AddConstraintCorrection(Constraint,AB,l_C,C)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: add constraint correction to coefficients.
C
      Implicit None
      Integer Constraint
      Integer AB
      Integer l_C
      Real*8  C(l_C)

      Character*27 SecNam
      Parameter (SecNam='LDF_AddConstraintCorrection')

      If (Constraint.eq.-1) Then ! no constraint
         Return
      Else If (Constraint.eq.0) Then ! charge constraint
         Call LDF_AddChargeConstraintCorrection(AB,l_C,C)
      Else ! unknown constraint
         Call WarningMessage(2,SecNam//': illegal constraint')
         Write(6,'(A,I10)') 'Constraint=',Constraint
         Call LDF_Quit(1)
      End If

      End
