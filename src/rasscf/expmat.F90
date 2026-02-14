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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Feb. 03, 2022, created this file.               *
! ****************************************************************
! Purpose: Return to the exponent of a skew-symmetric matrix
! Explanation:
! Calculate M=exp(-A)
! A is stored as a 1 by LenA vector, and M is a DimM by DimM matrix
! This subroutine assumes LenA=(DimM-1)*DimM/2
! How it works:
! First transform A into a skew-symmetric matrix form. Then call
! ExpMat_Inner to calculate exp(A)

subroutine ExpMat(M,A,DimM,LenA)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: DimM, LenA
real(kind=wp) :: M(DimM**2), A(LenA)
integer(kind=iwp) :: I, iIJ1, iIJ2, J, N
real(kind=wp), allocatable :: MatA(:)

call mma_allocate(MatA,DimM**2,Label='MatA')
N = DimM**2
call FZero(MatA,N)
do I=2,DimM
  do J=1,I-1
    iIJ2 = (I-2)*(I-1)/2+J
    iIJ1 = (I-1)*DimM+J
    MatA(iIJ1) = A(iIJ2)
    iIJ1 = (J-1)*DimM+I
    MatA(iIJ1) = -A(iIJ2)
  end do
end do

call ExpMat_Inner(M,MatA,DimM)
call mma_deallocate(MatA)

return

end subroutine ExpMat
