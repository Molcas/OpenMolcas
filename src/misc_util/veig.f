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
      SUBROUTINE VEIG(N,A,V)
C
C     THE DIAGONAL ELEMENTS OF THE LOWER TRIANGULAR MATRIX, A, STORED
C     IN UPPER PACKED STORAGE MODE ARE COPIED TO THE FIELD V.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),V(*)
C
      DO 10 I=1,N
       V(I)=A((I*(I+1))/2)
   10 CONTINUE
C
      RETURN
      END
