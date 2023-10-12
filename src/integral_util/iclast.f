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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************
      Integer Function iCLast(KWord,iChrct)
      Implicit None
      Character(LEN=*) KWord
      Integer iChrct

      Integer i

      iCLast = 0
      Do i = iChrct, 1, -1
         If (KWord(i:i).ne.' ') Then
            iCLast = i
            Return
         End If
      End Do

      End Function iCLast
