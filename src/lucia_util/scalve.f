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
      SUBROUTINE SCALVE(VECTOR,FACTOR,NDIM)
C
C CALCULATE SCALAR(FACTOR) TIMES VECTOR
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VECTOR(*)
C
      DO 100 I=1,NDIM
       VECTOR(I)=VECTOR(I)*FACTOR
  100 CONTINUE
C
      RETURN
      END
