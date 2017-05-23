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
      Real*8 FUNCTION DNORM2(N,X,INCX)
C
C     EUCLIDEAN NORM OF A VECTOR
C
      INTEGER    N, INCX
      REAL*8     X(*), DNRM2_
      External DNrm2_

      DNORM2 = DNRM2_(N,X,INCX)

      RETURN
      END
