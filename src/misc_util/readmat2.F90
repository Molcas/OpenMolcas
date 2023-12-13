!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on May. 19, 2022, created this file.               *
!*****************************************************************

subroutine ReadMat2(FileName,MatInfo,Matrix,NRow,NCol,LenName,LenInfo,Trans)
! This subroutine is to replace ReadMat in the long run.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NRow, NCol, LenName, LenInfo
character(len=LenName), intent(in) :: FileName
character(len=LenInfo), intent(out) :: MatInfo
real(kind=wp), intent(out) :: Matrix(NRow*NCol)
character, intent(in) :: Trans
integer(kind=iwp) :: ICol, iOff, IRow, LU
integer(kind=iwp), external :: IsFreeUnit

if (LenName > 0) then
  LU = IsFreeUnit(100)
  call Molcas_Open(LU,FileName)
else
  LU = u6
end if
if (Trans == 'T') then
  do IRow=1,NRow
    iOff = (IRow-1)*nCol
    read(LU,*) (Matrix(iOff+ICol),ICol=1,NCol)
  end do
else
  do ICol=1,NCol
    read(LU,*) (Matrix((iRow-1)*nCol+iCol),IRow=1,NRow)
  end do
end if
read(LU,*) MatInfo
if (LenName > 0) close(LU)

return

end subroutine ReadMat2
