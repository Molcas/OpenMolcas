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
      REAL*8 FUNCTION INPROD_MCLR(A,B,NDIM)
C      CALCULATE SCALAR PRODUCT BETWEEN TO VECTORS A,B
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*)
      INPROD_MCLR=0.0D0
      DO 100 I=1,NDIM
       INPROD_MCLR=INPROD_MCLR+A(I)*B(I)
  100 CONTINUE
C
      RETURN
      END
