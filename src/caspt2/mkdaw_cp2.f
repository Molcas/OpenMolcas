************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKDAW_CP2(NLEV,NVERT,IDRT,IDOWN,IDAW,LTV)
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer NLEV, NVERT
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4),IDRT(NVERT,5)
      DIMENSION LTV(-1:NLEV)
      PARAMETER (LTAB=1)

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

#ifdef _DEBUGPRINT_
C CHECK PRINTS:
      WRITE(6,*)
      WRITE(6,*)' LEVEL-TO-VERTEX TABLE:'
      DO LEV=0,NLEV
        WRITE(6,'(1X,I4,5X,I4,A4,I4)') LEV,LTV(LEV),' -- ',LTV(LEV-1)-1
      END DO
#endif
      END SUBROUTINE MKDAW_CP2
