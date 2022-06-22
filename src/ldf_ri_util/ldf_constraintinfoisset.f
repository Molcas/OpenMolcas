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
      Logical Function LDF_ConstraintInfoIsSet(Constraint)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: return .True. if constraint info is set.
C
      Implicit None
      Integer Constraint

      Logical LDF_ChargeConstraintInfoIsSet

      If (Constraint.eq.-1) Then ! no constraint
         LDF_ConstraintInfoIsSet=.True.
      Else If (Constraint.eq.0) Then ! charge constraint
         LDF_ConstraintInfoIsSet=LDF_ChargeConstraintInfoIsSet()
      Else ! unknown constraint
         Call WarningMessage(2,
     &                    'LDF_ConstraintInfoIsSet: unknown constraint')
         Call LDF_Quit(1)
         LDF_ConstraintInfoIsSet=.False. ! avoid compiler warnings
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Logical Function LDF_ChargeConstraintInfoIsSet()
      Implicit None
#include "ldf_charge_constraint_info.fh"
      LDF_ChargeConstraintInfoIsSet=ChargeConstraintSet
      End
