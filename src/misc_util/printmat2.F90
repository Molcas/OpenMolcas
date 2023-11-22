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

subroutine PrintMat2(FileName,MatInfo,Matrix,NRow,NCol,LenName,LenInfo,Trans)
! This subroutine is to replace PrintMat in the long run.
! Matrix is now a nRow*nCol array.
! Note that the column index is the fast running index in Fortran,
! so when TRANS='T', it prints the matrix by proceeding with the
! fast-running index.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NRow, NCol, LenName, LenInfo
character(len=LenName), intent(in) :: FileName
character(len=LenInfo), intent(in) :: MatInfo
real(kind=wp), intent(in) :: Matrix(NRow*NCol)
character, intent(in) :: Trans
character(len=80) :: PrtFmt
integer(kind=iwp) :: ICol, IRow, iOff, LU
integer(kind=iwp), external :: IsFreeUnit

if (LenName > 0) then
  LU = IsFreeUnit(100)
  call Molcas_Open(LU,FileName)
else
  LU = u6
end if
if (Trans == 'T') then
  write(PrtFmt,'(A,I5,A)') '(',NCol,'(ES24.14E4,1X))'
  do IRow=1,NRow
    iOff = (IRow-1)*NCol
    write(LU,PrtFmt) (Matrix(iOff+ICol),ICol=1,NCol)
  end do
else
  write(PrtFmt,'(A,I5,A)') '(',NRow,'(ES24.14E4,1X))'
  do ICol=1,NCol
    write(LU,PrtFmt) (Matrix((iRow-1)*NCol+iCol),IRow=1,NRow)
  end do
end if
write(LU,*) MatInfo
if (LenName > 0) close(LU)

return

end subroutine PrintMat2
