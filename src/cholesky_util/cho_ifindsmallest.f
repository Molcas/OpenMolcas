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
* Copyright (C) 2007, Thomas Bondo Pedersen                            *
************************************************************************
      Integer Function Cho_iFindSmallest(iVec,n)
C
C     Thomas Bondo Pedersen, June 2007.
C
C     Return index of smallest element in array iVec.
C
      Implicit None
      Integer n
      Integer iVec(n)

      Integer i, imin

      If (n .lt. 1) Then
         imin = 0
      Else
         imin = n
         Do i = n-1,1,-1
            If (iVec(i) .lt. iVec(imin)) Then
               imin = i
            End If
         End Do
      End If

      Cho_iFindSmallest = imin

      End
