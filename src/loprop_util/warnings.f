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
      Subroutine Warnings(iCode,Warning,iLength)
      Implicit Real*8 (A-H,O-Z)
      Character*(*) Warning
*
* Binary codes:
*  1: Multiple minima found
*  2: Minima not within range
*  3: Ran out of iterations
*
      If (iLength .lt. 25) Then
         Write(6,*)
     &         'Length of warning string must be at least 25 characters'
         Call Abend()
      End If
*
      Do i = 1,iLength
         Warning(i:i) = ' '
      End Do
*
      If (iCode .eq. 1) Then
         Warning = 'Multiple minima found'
      Else If (iCode .eq. 2) Then
         Warning = 'Minima not within range'
      Else If (iCode .eq. 3) Then
         Warning = 'Ran out of iterations'
      Else If (iCode .eq. 4) Then
         Warning = 'No minima found'
      End If
*
      Return
      End
