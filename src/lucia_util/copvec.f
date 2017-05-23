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
      SUBROUTINE COPVEC(FROM,TO,NDIM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/COPVECST/XNCALL_COPVEC, XNMOVE_COPVEC
      DIMENSION FROM(*),TO(*)
C
      XNCALL_COPVEC = XNCALL_COPVEC + 1.0D0
      XNMOVE_COPVEC = XNMOVE_COPVEC + DBLE(NDIM)
      DO 100 I=1,NDIM
       TO(I)=FROM(I)
  100 CONTINUE
C
      RETURN
      END
