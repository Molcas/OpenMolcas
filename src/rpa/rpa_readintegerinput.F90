!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine RPA_ReadIntegerInput(Key,nInp,Lu,iVal,n)
      Implicit None
      Character*4 Key
      Integer nInp
      Integer Lu
      Integer n
      Integer iVal(n)
      Character*180 Line

      Character*180 Get_Ln
      External Get_Ln

      If (n.ge.nInp) Then
         Line=Get_Ln(Lu)
         Call Get_I(1,iVal,nInp)
      Else
         ! insufficent memory for reading (fix in calling routine)
         Call RPA_Warn(3,'Integer read problem for keyword '//Key)
      End If

      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_character(Line)
#endif
      End
