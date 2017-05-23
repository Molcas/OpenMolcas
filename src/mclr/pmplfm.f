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
      SUBROUTINE PMPLFM(AP,B,NDIM)
*
* Add lower half of a full matrix to a matrix packed
* in lower triangular form ( packed matrix stored columnwise )
*
      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION AP(*),B(*)
*
      IBSP = 1
      IBSF = 1
      DO 100 ICOL = 1, NDIM
        NELMNT = NDIM - ICOL + 1
        CALL VECSUM(AP(IBSP),AP(IBSP),B(IBSF),1.0D0,1.0D0,NELMNT)
        IBSP = IBSP + NELMNT
        IBSF = IBSF + NDIM
  100 CONTINUE
*
      RETURN
      END
