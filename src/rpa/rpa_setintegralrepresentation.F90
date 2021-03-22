!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine RPA_SetIntegralRepresentation()
      Implicit None
#include "rpa_config.fh"
      Call DecideOnCholesky(doCD)
      Call DecideOnDF(doDF)
      Call DecideOnLocalDF(doLDF)
      If (doLDF) Then
         doCD=.false.
         doDF=.false.
      Else If (doDF) Then
         doCD=.false.
         doLDF=.false.
      Else If (doCD) Then
         doDF=.false.
         doLDF=.false.
      End If
      End
