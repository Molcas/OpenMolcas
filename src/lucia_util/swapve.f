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
      SUBROUTINE SWAPVE(VEC1,VEC2,NDIM)
C
C      SWAP ELEMENTS OF VECTORS VEC1 AND VEC2
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VEC1(*) ,VEC2(*)
      DO 100 I=1,NDIM
       BUF=VEC1(I)
       VEC1(I)=VEC2(I)
       VEC2(I)=BUF
  100 CONTINUE
C
      RETURN
      END
