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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine Copy2DMat(A,B,NRow,NCol)

implicit none
integer NRow, NCol, IRow, ICol
real*8, dimension(NRow,NCol) :: A, B

do ICol=1,NCol
  do IRow=1,NRow
    A(IRow,ICol) = B(IRow,ICol)
  end do
end do

end subroutine Copy2DMat
