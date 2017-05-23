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
      SUBROUTINE MKDAW_CP2(IDRT,IDOWN,IDAW,LTV)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "pt2_guga.fh"
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
      DO IC=0,3
        IDAW(NVERT,IC)=0
      END DO
      IDAW(NVERT,4)=1
      DO IV=NVERT-1,1,-1
        ISUM=0
        DO IC=0,3
          IDAW(IV,IC)=0
          IDWN=IDOWN(IV,IC)
          IF(IDWN.EQ.0) GOTO 20
          IDAW(IV,IC)=ISUM
          ISUM=ISUM+IDAW(IDWN,4)
  20      CONTINUE
        END DO
        IDAW(IV,4)=ISUM
      END DO

#ifdef _DEBUG_
C CHECK PRINTS:
      WRITE(6,*)
      WRITE(6,*)' LEVEL-TO-VERTEX TABLE:'
      DO LEV=0,NLEV
        WRITE(6,'(1X,I4,5X,I4,A4,I4)') LEV,LTV(LEV),' -- ',LTV(LEV-1)-1
      END DO
      WRITE(6,*)
      WRITE(6,*)' DIRECT ARC WEIGHTS:'
      DO IV=1,NVERT
        WRITE(6,'(1X,I4,5X,5(1X,I6))') IV,(IDAW(IV,IC),IC=0,4)
      END DO
      WRITE(6,*)
#endif
      RETURN
      END
