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

subroutine Cho_RdQCol_Indx(xInt,IDCol,nRow,nCol,Lunit)
!
! Purpose: read indexed qualified columns from WA file with unit
!          Lunit. (WA=word-addressable)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCol, IDCol(nCol), nRow, Lunit
real(kind=wp), intent(out) :: xInt(nRow,nCol)
integer(kind=iwp) :: iAdr, iCol, iOpt, lTot

if ((nRow < 1) .or. (nCol < 1)) return

do iCol=1,nCol
  iOpt = 2
  lTot = nRow
  iAdr = nRow*(IDCol(iCol)-1)
  call dDAFile(Lunit,iOpt,xInt(:,iCol),lTot,iAdr)
end do

end subroutine Cho_RdQCol_Indx
