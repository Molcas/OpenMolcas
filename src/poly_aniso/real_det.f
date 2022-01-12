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
      REAL*8 FUNCTION FindDetR(matrix, n)
c     Function to find the determinant of a square matrix
c
c     Description: The Subroutine is based on two key points:
c
c     1]  A determinant is unaltered when row operations are performed:
c         Hence, using this principle, row operations (column operations
c         would work as well) are used to convert the matrix the matrix
c         into upper triangular form
c     2]  The determinant of a triangular matrix is obtained by finding
c         the product of the diagonal elements
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: N
      REAL (kind=8)        :: matrix(N,N)
      ! local variables:
      REAL (kind=8)        :: m, temp, MINIMAL_REAL
      INTEGER              :: i, j, k, l
      LOGICAL              :: DetExists
      DetExists = .TRUE.
      MINIMAL_REAL = TINY(0.d0)
      l = 1
      temp=0
      ! Convert to upper triangular form
      DO k = 1, N-1
         IF (ABS(matrix(k,k))<MINIMAL_REAL) THEN
            DetExists = .FALSE.
            DO i = k+1, N
               IF (ABS(matrix(i,k))>MINIMAL_REAL) THEN
                  DO j = 1, N
                           temp  =  matrix(i,j)
                     matrix(i,j) =  matrix(k,j)
                     matrix(k,j) =  temp
                  END DO
                  DetExists = .TRUE.
               END IF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
               FindDetR = 0
               RETURN
            END IF
         END IF
         DO j = k+1, N
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, N
               matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
         END DO
      END DO ! k
      ! Evaluate determinant by finding product of diagonal elements
      FindDetR = l
      DO i=1, N
         FindDetR = FindDetR * matrix(i,i)
      END DO
      RETURN
      END FUNCTION FindDetR
