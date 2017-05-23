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
      SubRoutine ChoMP2_Col_cp(X,nRow,nCol,S,nSRow,iSRow)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: copy subblock of matrix.
C
      Implicit None
      Integer nRow, nCol, nSRow
      Integer iSRow(nSRow)
      Real*8  X(nRow,nCol), S(nSRow,nCol)

      Integer iCol, iSR, iR

      Do iCol = 1,nCol
         Do iSR = 1,nSRow
            iR = iSRow(iSR)
            S(iSR,iCol) = X(iR,iCol)
         End Do
      End Do

      End
