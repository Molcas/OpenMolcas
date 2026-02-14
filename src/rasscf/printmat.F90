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

subroutine PrintMat(FileName,MatInfo,Matrix,NRow,NCol,LenName,LenInfo,Trans)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NRow, NCol, LenName, LenInfo
character(len=LenName) :: FileName
character(len=LenInfo) :: MatInfo
real(kind=wp) :: Matrix(NRow,NCol)
integer(kind=iwp) :: ICol, IRow, LU
character(len=80) :: PrtFmt
character(len=1) :: Trans
integer(kind=iwp), external :: IsFreeUnit

if (LenName > 0) then
  LU = IsFreeUnit(100)
  call Molcas_Open(LU,FileName)
else
  LU = u6
end if

if (Trans == 'N') then
  write(PrtFmt,'(A,I5,A)') '(',NCol,'(ES24.14E4,1X))'
  do IRow=1,NRow
    write(LU,PrtFmt) (Matrix(IRow,ICol),ICol=1,NCol)
  end do
else
  write(PrtFmt,'(A,I5,A)') '(',NRow,'(ES24.14E4,1X))'
  do ICol=1,NCol
    write(LU,PrtFmt) (Matrix(IRow,ICol),IRow=1,NRow)
  end do
end if
write(LU,*) MatInfo
if (LenName > 0) close(LU)

end subroutine PrintMat
