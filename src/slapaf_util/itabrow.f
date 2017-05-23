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
      Integer Function iTabRow(i)
      Implicit Integer (a-z)
*
      iTabRow=0
      iTabRow=1
      If (i.gt. 0 .and. i.le. 2) Then
         iTabRow=1
      Else If (i.gt. 2 .and. i.le.10) Then
         iTabRow=2
      Else If (i.gt.10 .and. i.le.18) Then
         iTabRow=3
      Else If (i.gt.18 .and. i.le.36) Then
         iTabRow=4
      Else If (i.gt.36 .and. i.le.54) Then
         iTabRow=5
      Else If (i.gt.54 .and. i.le.86) Then
         iTabRow=6
      Else If (i.gt.86) Then
         iTabRow=7
      End If
*
      Return
      End
