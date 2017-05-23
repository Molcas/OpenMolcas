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
      SUBROUTINE ZEROMA(W,I1,I2)
      implicit none
      integer I1,I2, I
      REAL*8 ZERO,W
      PARAMETER (ZERO=0.D0)
      DIMENSION W(*)
      IF(I2.LT.I1)RETURN
      DO I=I1,I2
         W(I)=ZERO
      ENDDO
      RETURN
      END
