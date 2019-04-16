************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C
C     Computes eigenvalues and optionally (right) eigenvectors of a
C     general square matrix
C
      SUBROUTINE XEIGEN (NVEC,NA,N,A,EVR,EVI,VECS,IERR)
      INTEGER NVEC,NA,N,IERR
      REAL*8 A(NA,N),EVR(N),EVI(N),VECS(NA,N)
*
#include "stdalloc.fh"
      CHARACTER JL,JR
      INTEGER NW
      REAL*8 TMP(1)
      REAL*8, DIMENSION(:), ALLOCATABLE :: WRK
      JL='N'
      JR='N'
      IF (NVEC.NE.0) JR='V'
      IERR=0
      CALL DGEEV_(JL,JR,N,A,NA,EVR,EVI,VECS,NA,VECS,NA,TMP,-1,IERR)
      NW=INT(TMP(1))
      CALL mma_allocate(WRK,NW)
      CALL DGEEV_(JL,JR,N,A,NA,EVR,EVI,VECS,NA,VECS,NA,WRK,NW,IERR)
      CALL mma_deallocate(WRK)
      END
