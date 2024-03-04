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
      SUBROUTINE mkRAS(SGS)
      use Struct, only: SGStruct
      use UnixInfo, only: ProgName
      IMPLICIT None
      Type(SGStruct) SGS

      If (ProgName(1:5).eq.'rassi') Then
         Call rmvert(SGS)
      Else
         Call RESTR(SGS)
      End If

      END SUBROUTINE mkRAS
