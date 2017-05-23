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
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM,ISING)
C FIND INVERSE OF MATRIX A
C INPUT :
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
*
* ISING = 0 => No convergence problems
*       = 1  => Convergence problems
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      ITEST=0
      IF(NDIM.EQ.1)THEN
        IF(A(1,1) .NE. 0.0D0 ) THEN
           A(1,1) = 1.0D0/A(1,1)
        ELSE
           ITEST = 1
        END IF
      ELSE
        DETERM=0.0D0
        EPSIL=0.0D0
        CALL BNDINV(        A,        B,     NDIM,   DETERM,    EPSIL,
     &                  ITEST,   MATDIM)
      END IF
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
*
      IF(ITEST.NE.0) THEN
        ISING = 1
      ELSE
        ISING = 0
      END IF
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
