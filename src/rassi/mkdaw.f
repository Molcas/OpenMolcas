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
      SUBROUTINE MKDAW_RASSI(NLEV,NVERT,IDRT,IDOWN,IDAW,LTV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4),IDRT(NVERT,5)
      DIMENSION LTV(-1:NLEV)
      PARAMETER (LTAB=1)

C Purpose: Construct Direct Arc Weights.

C SET UP A LEVEL-TO-VERTEX TABLE, LTV, AND IDENTIFY MIDVERTICES:
      DO LEV=-1,NLEV
        LTV(LEV)=0
      END DO
      DO IV=1,NVERT
        LEV=IDRT(IV,LTAB)
        LTV(LEV)=LTV(LEV)+1
      END DO
      DO LEV=NLEV,0,-1
        LTV(LEV-1)=LTV(LEV-1)+LTV(LEV)
      END DO
      DO LEV=-1,NLEV-1
        LTV(LEV)=1+LTV(LEV+1)
      END DO

      CALL MKDAW(NVERT,IDOWN,IDAW)

      END
