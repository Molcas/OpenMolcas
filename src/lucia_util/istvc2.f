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
      SUBROUTINE ISTVC2(IVEC,IBASE,IFACT,NDIM)
C
C IVEC(I) = IBASE + IFACT * I
C
C     DIMENSION IVEC(1   )
      DIMENSION IVEC(NDIM)
C
      DO 100 I = 1,NDIM
        IVEC(I) = IBASE + IFACT*I
  100 CONTINUE
C
      RETURN
      END
