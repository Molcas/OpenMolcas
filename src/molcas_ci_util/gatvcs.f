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

      SUBROUTINE GATVCS(VECO,VECI,INDEX,NDIM)
* Gather vector allowing for sign change
*
* VECO(I) = VECI(INDEX(I))
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION VECI(*),VECO(*),INDEX(*)
      INTRINSIC SIGN
*
      DO 100 I = 1, NDIM
      VECO(I) = VECI(ABS(INDEX(I)))*DBLE(SIGN(1,INDEX(I)))
  100 CONTINUE
*
      RETURN
      END
      SUBROUTINE SCAVCS(VECO,VECI,INDEX,NDIM)
*
* Scatter vector with sign change
*
* vecO(abs(index(i))) = veci(i)*sign(index(i))
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION VECI(*),VECO(*),INDEX(*)
      INTRINSIC SIGN
C
      DO 100 I = 1, NDIM
      VECO(ABS(INDEX(I))) = VECI(I)*DBLE(SIGN(1,INDEX(I)))
  100 CONTINUE
C
      RETURN
      END
