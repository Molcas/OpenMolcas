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

subroutine XEIGEN(NVEC,NA,N,A,EVR,EVI,VECS,IERR)
! Computes eigenvalues and optionally (right) eigenvectors of a
! general square matrix

integer NVEC, NA, N, IERR
real*8 A(NA,N), EVR(N), EVI(N), VECS(NA,N)

#include "stdalloc.fh"
character JL, JR
integer NW
real*8 TMP(1)
real*8, dimension(:), allocatable :: WRK

JL = 'N'
JR = 'N'
if (NVEC /= 0) JR = 'V'
IERR = 0
call DGEEV_(JL,JR,N,A,NA,EVR,EVI,VECS,NA,VECS,NA,TMP,-1,IERR)
NW = int(TMP(1))
call mma_allocate(WRK,NW)
call DGEEV_(JL,JR,N,A,NA,EVR,EVI,VECS,NA,VECS,NA,WRK,NW,IERR)
call mma_deallocate(WRK)

end subroutine XEIGEN
