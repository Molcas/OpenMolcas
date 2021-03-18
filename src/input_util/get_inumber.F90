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
! Copyright (C) 2004, Giovanni Ghigo                                   *
!***********************************************************************
!***********************************************************************
! This routine checks if -string- is an integer. If it is a number     *
! NaN=.False. and -iNumber- contains the number, otherwise, NaN=.True. *
!----------------------------------------------------------------------*
! Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
! Author:   Giovanni Ghigo                                             *
!***********************************************************************
      Subroutine Get_iNumber(string,iNumber,iErr)
      Implicit Integer (i-n)
      Character*(*)  string
      Character*11  Chars
      Data Chars /' 1234567890' /
      Logical NaN
      iErr=0
      iNumber=0
      NaN=.True.
      Do i=1,len(string)
      NaN=.True.
        Do j=1,11
          If (string(i:i).EQ.Chars(j:j)) NaN=.False.
        EndDo
      If (NaN) GoTo 10
      End Do
10    If (.NOT.NaN) Read(string,*) iNumber
      If (NaN) iErr=1
      Return
      End
