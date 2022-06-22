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
      Subroutine LDF_SetConstraint(Constraint)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: set info for constrained LDF.
C              To unset, call LDF_UnsetConstraint()
C
      Implicit None
      Integer Constraint

      Character*17 SecNam
      Parameter (SecNam='LDF_SetConstraint')

      If (Constraint.eq.-1) Then ! no constraint
         Return
      Else If (Constraint.eq.0) Then ! charge constraint
         Call LDF_SetChargeConstraint()
      Else ! unknown constraint
         Call WarningMessage(2,SecNam//': illegal constraint')
         Write(6,'(A,I10)') 'Constraint=',Constraint
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetConstraint(Constraint)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: unset info for constrained LDF.
C
      Implicit None
      Integer Constraint

      Character*17 SecNam
      Parameter (SecNam='LDF_SetConstraint')

      If (Constraint.eq.-1) Then ! no constraint
         Return
      Else If (Constraint.eq.0) Then ! charge constraint
         Call LDF_UnsetChargeConstraint()
      Else ! unknown constraint
         Call WarningMessage(2,SecNam//': illegal constraint')
         Write(6,'(A,I10)') 'Constraint=',Constraint
         Call LDF_Quit(1)
      End If

      End
