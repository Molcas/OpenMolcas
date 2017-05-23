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
      Subroutine Set_Basis_Mode(Label)
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
      Character*(*) Label
      Character*7 Lbl
*
      Atomic=.False.
      kCnttp=0
      Lbl=Label(1:7)
      Call UpCase(Lbl)
      If (Lbl.eq.'VALENCE') Then
         Basis_Mode = Valence_Mode
      Else If (Lbl.eq.'AUXILIA') Then
         Basis_Mode = Auxiliary_Mode
      Else If (Lbl.eq.'FRAGMEN') Then
         Basis_Mode = Fragment_Mode
      Else If (Lbl.eq.'WITHAUX') Then
         Basis_Mode = With_Auxiliary_Mode
      Else If (Lbl.eq.'WITHFRA') Then
         Basis_Mode = With_Fragment_Mode
      Else If (Lbl.eq.'ALL    ') Then
         Basis_Mode = All_Mode
      Else
         Call WarningMessage(2,'Set_Basis_Mode: illegal mode,'//
     &               'Label='//Lbl)
         Call Abend()
      End If
*
      Return
      End
