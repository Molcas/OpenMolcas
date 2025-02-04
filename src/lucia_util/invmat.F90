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

subroutine INVMAT(A,B,MATDIM,NDIM,ISING)
! FIND INVERSE OF MATRIX A
! INPUT :
!        A : MATRIX TO BE INVERTED
!        B : SCRATCH ARRAY
!        MATDIM : PHYSICAL DIMENSION OF MATRICES
!        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
!
! OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
! WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
!
! ISING = 0 => No convergence problems
!       = 1 => Convergence problems

use Constants, only: Zero, One
use Definitions, only: u6

implicit real*8(A-H,O-Z)
dimension A(MATDIM,MATDIM), B(MATDIM,MATDIM)

ITEST = 0
if (NDIM == 1) then
  if (A(1,1) /= Zero) then
    A(1,1) = One/A(1,1)
  else
    ITEST = 1
  end if
else
  DETERM = Zero
  EPSIL = Zero
  call BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
end if

if (ITEST /= 0) then
  write(u6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
  ISING = 1
else
  ISING = 0
end if

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' INVERTED MATRIX'
  call WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
end if

end subroutine INVMAT
