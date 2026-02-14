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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine ReadMat(FileName,MatInfo,Matrix,NRow,NCol,LenName,LenInfo,Trans)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NRow, NCol, LenName, LenInfo
character(len=LenName) :: FileName
character(len=LenInfo) :: MatInfo
real(kind=wp) :: Matrix(NRow,NCol)
character :: Trans
integer(kind=iwp) :: ICol, IRow, LU
integer(kind=iwp), external :: IsFreeUnit

if (LenName > 0) then
  LU = IsFreeUnit(100)
  call Molcas_Open(LU,FileName)
else
  LU = u6
end if
if (Trans == 'N') then
  do IRow=1,NRow
    read(LU,*) (Matrix(IRow,ICol),ICol=1,NCol)
  end do
else
  do ICol=1,NCol
    read(LU,*) (Matrix(IRow,ICol),IRow=1,NRow)
  end do
end if
read(LU,*) MatInfo
if (LenName > 0) close(LU)

end subroutine ReadMat
