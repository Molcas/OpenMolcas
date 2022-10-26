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

integer NRow, NCol, LenName
character(Len=LenName) :: FileName
character(Len=LenInfo) :: MatInfo
character(Len=1) :: Trans
character(Len=80) :: PrtFmt
real*8, dimension(NRow*NCol) :: Matrix
integer LU, IsFreeUnit, IRow, ICol, iOff
external IsFreeUnit

if (LenName > 0) then
  LU = 100
  LU = IsFreeUnit(LU)
  call Molcas_Open(LU,FileName)
else
  LU = 6
end if
if (Trans == 'T') then
  write(PrtFmt,'(A1,I5,A14)') '(',NCol,'(E24.14E4,1X))'
  do IRow=1,NRow
    iOff = (IRow-1)*nCol
    write(LU,PrtFmt) (Matrix(iOff+ICol),ICol=1,NCol)
  end do
else
  write(PrtFmt,'(A1,I5,A14)') '(',NRow,'(E24.14E4,1X))'
  do ICol=1,NCol
    write(LU,PrtFmt) (Matrix((iRow-1)*nCol+iCol),IRow=1,NRow)
  end do
end if
write(LU,*) MatInfo
if (LenName > 0) then
  close(LU)
end if

return

end subroutine PrintMat2
