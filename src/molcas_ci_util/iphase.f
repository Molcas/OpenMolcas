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
C
C
      FUNCTION IPHASE(IDRT,IUP,IWALK)
C
C     PURPOSE: THE SYMMETRIC GROUP APPROACH AND THE UNITARY GROUP
C              APPROACH DIFFER IN THE PHASE CONVETION. FIND THE
C              PHASE FACTOR RELATING THE CSFS IN EITHER BASIS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "gugx.fh"
C
      DIMENSION IDRT(NVERT,5),IUP(NVERT,0:3),IWALK(NLEV)
C
C     FIND THE MIDVERTEX AND THE COMBINED WALK SYMMETRY
C
      IPHASE=1
      IVERT=NVERT
      DO 100 LEV=1,NLEV
        ICASE=IWALK(LEV)
        IVERT=IUP(IVERT,ICASE)
        IF( ICASE.EQ.2 .OR. ICASE.EQ.3 ) THEN
          ISIGN=(-1)**IDRT(IVERT,4)
          IPHASE=IPHASE*ISIGN
        ENDIF
100   CONTINUE
C
C     EXIT
C
      RETURN
      END
