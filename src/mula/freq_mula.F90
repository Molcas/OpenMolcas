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
! Copyright (C) 1994, Niclas Forsberg                                  *
!***********************************************************************

subroutine Freq_mula(Hess,G,V,Lambda,B,qMat,NumInt,NumOfAt)
!  Purpose:
!    Find eigenvalues and eigenvectors of FG matrix.
!    The eigenvalues are stored in the array Lambda and the eigen-
!    vectors are stored as the columns of the matrix V.
!    The eigenvectors are then used to calculate the cartesian
!    displacements, which are stored in matrix X.
!
!  Input:
!    Hess     : Real two dimensional array -  contains the force constants expressed in internal
!    G        : Real two dimensional array.
!    B        : Real two dimensional array.
!
!  Output:
!    V        : Real two dimensional array - contains the eigenvectors of F*G as columns.
!    Lambda   : Real array - contains the eigenvalues of F*G.
!    qMat     : Real array - contains the cartesian displacements of the atoms.
!
!  Local:
!    U,Tmp    : Real two dimensional arrays.
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NumInt, NumOfAt
real(kind=wp), intent(in) :: Hess(NumInt,NumInt), G(NumInt,NumInt), B(3*NumOfAt,NumInt)
real(kind=wp), intent(out) :: V(NumInt,NumInt), Lambda(NumInt), qMat(3*NumOfAt,NumInt)
real(kind=wp) :: Det
real(kind=wp), allocatable :: Temp(:,:), U(:,:)

! Solve secular equation.
call SolveSecEq(Hess,NumInt,V,G,Lambda)

! get memory space for U.
call mma_allocate(U,NumInt,NumInt,label='U')

! Copy matrix containing eigenvectors to U, because this matrix
! will be destroyed when subroutine Dool_MULA is called.
U(:,:) = V

! Calculate the cartesian diplacements, i.e. solve
!     qMat = B * ( B(T) * B )^(-1) * V.

! get memory for matrix Temp.
call mma_allocate(Temp,NumInt,NumInt,label='Temp')

call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,B,3*NumOfAt,B,3*NumOfAt,Zero,Temp,NumInt)
call Dool_MULA(Temp,NumInt,NumInt,U,NumInt,NumInt,Det)
call DGEMM_('N','N',3*NumOfAt,NumInt,NumInt,One,B,3*NumOfAt,U,NumInt,Zero,qMat,3*NumOfAt)

! free memory space of Temp and U.
call mma_deallocate(U)
call mma_deallocate(Temp)

end subroutine Freq_mula
