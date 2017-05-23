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
      SUBROUTINE PRDRT_m(NVERT,DRT,DOWN)
C
C     PURPOSE: PRINT THE DRT TABLE
C
      IMPLICIT INTEGER (A-Z)
#include "output_ras.fh"
      Parameter (ROUTINE='PRDRT   ')
C
      DIMENSION DRT(NVERT,5),DOWN(NVERT,0:3)
C
      Write(LF,*)
      Write(LF,*)' VERT      L  N    A  B  C      CHAINING INDICES.'
      DO V=1,NVERT
        Write(LF,'(1X,I4,5X,2I3,2X,3I3,5X,4I4)')
     &              V,(DRT(V,I),I=1,5),(DOWN(V,S),S=0,3)
      END DO
      Write(LF,*)
C
C     EXIT
C
      RETURN
      END
