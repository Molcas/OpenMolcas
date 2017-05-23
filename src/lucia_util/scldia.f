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
      SUBROUTINE SCLDIA(A,FACTOR,NDIM,IPACK)
*
* scale diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*                assumed packed columnwise !!!!
      IMPLICIT REAL*8 (A-H,O-Z)
*
      DIMENSION A(*)
*
      IF( IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          II = (I-1)*NDIM + I
          A(II) = A(II) * FACTOR
  100   CONTINUE
      ELSE
        II = 1
        DO 200 I = 1, NDIM
          A(II) = A(II) * FACTOR
          II = II + NDIM - I + 1
  200   CONTINUE
      END IF
*
      RETURN
      END
