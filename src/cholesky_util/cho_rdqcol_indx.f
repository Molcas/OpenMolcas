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
      SubRoutine Cho_RdQCol_Indx(xInt,IDCol,nRow,nCol,Lunit)
C
C     Purpose: read indexed qualified columns from WA file with unit
C              Lunit. (WA=word-addressable)
C
      Implicit None
      Integer nRow, nCol
      Real*8  xInt(nRow,nCol)
      Integer IDCol(nCol)
      Integer Lunit

      Integer iCol, iOpt, iAdr, lTot

      If (nRow.lt.1 .or. nCol.lt.1) Return

      Do iCol = 1,nCol
         iOpt = 2
         lTot = nRow
         iAdr = nRow*(IDCol(iCol)-1)
         Call dDAFile(Lunit,iOpt,xInt(1,iCol),lTot,iAdr)
      End Do

      End
