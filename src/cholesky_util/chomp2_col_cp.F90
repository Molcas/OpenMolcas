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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Col_cp(X,nRow,nCol,S,nSRow,iSRow)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: copy subblock of matrix.

implicit none
integer nRow, nCol, nSRow
integer iSRow(nSRow)
real*8 X(nRow,nCol), S(nSRow,nCol)
integer iCol, iSR, iR

do iCol=1,nCol
  do iSR=1,nSRow
    iR = iSRow(iSR)
    S(iSR,iCol) = X(iR,iCol)
  end do
end do

end subroutine ChoMP2_Col_cp
