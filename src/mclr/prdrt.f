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
      SUBROUTINE PRDRT_MCLR(NVERT,DRT,DOWN)
      IMPLICIT INTEGER (A-Z)
      DIMENSION DRT(NVERT,5),DOWN(NVERT,0:3)
      WRITE(6,*)
      CALL XFLUSH(6)
      WRITE(6,*)' VERT      L  N    A  B  C      CHAINING INDICES.'
      CALL XFLUSH(6)
      DO 10 V=1,NVERT
        WRITE(6,1000) V,(DRT(V,I),I=1,5),(DOWN(V,S),S=0,3)
      CALL XFLUSH(6)
10    CONTINUE
      WRITE(6,*)
      CALL XFLUSH(6)
      RETURN
1000  FORMAT(1X,I4,5X,2I3,2X,3I3,5X,4I4)
      END
