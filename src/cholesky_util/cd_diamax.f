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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine CD_DiaMax(Diag,nDim,iPivot,iQual,nQual,DiaMin)
C
C     Thomas Bondo Pedersen, October 2004.
C
C     Purpose: find nQual largest elements in Diag and leave pointers to
C              them in iQual. Only elements larger than DiaMin are
C              returned (thus, nQual may be reduced here).
C
      Implicit None
      Integer nDim, nQual
      Real*8  Diag(nDim)
      Real*8  DiaMin
      Integer iPivot(nDim), iQual(nQual)

      Integer i, j, iTmp, iMax

      Do i = 1,nDim
         iPivot(i) = i
      End Do

      Do j = 1,nQual
         Do i = nDim,j+1,-1
            If (Diag(iPivot(i)) .gt. Diag(iPivot(i-1))) Then
               iTmp = iPivot(i-1)
               iPivot(i-1) = iPivot(i)
               iPivot(i)   = iTmp
            End If
         End Do
      End Do

      Call iZero(iQual,nQual)
      iMax  = nQual
      i     = 0
      nQual = 0
      Do While (i .lt. iMax)
         i = i + 1
         If (Diag(iPivot(i)) .ge. DiaMin) Then
            nQual = nQual + 1
            iQual(nQual) = iPivot(i)
         Else
            i = iMax + 1
         End If
      End Do

      End
