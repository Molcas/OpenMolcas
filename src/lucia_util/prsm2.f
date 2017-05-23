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
      SUBROUTINE PRSM2(A,NDIM)
C
C PRINT LOWER TRIANGULAR MATRIX PACKED IN COLUMN WISE FASHION
C
      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION A(*)
C
      DO 100 I=1,NDIM
        WRITE(6,1010) I,
     &  (A((J-1)*NDIM-J*(J-1)/2+I),J=1,I)
  100 CONTINUE
      RETURN
 1010 FORMAT(1H0,2X,I3,5(E14.7),/,(1H ,5X,5(E14.7)))
      END
