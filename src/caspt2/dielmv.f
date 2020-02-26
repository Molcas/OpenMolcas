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
      SUBROUTINE DIELMV(ICASE,JCASE,NUP,NDWN,EMU)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ICASE(*),JCASE(*)
      DIMENSION EMU(NUP,NDWN)
#include "pt2_guga.fh"

      DO I=1,NUP
        II=NIPWLK*(I-1)
        SUM=0.0D0
        DO LV1=MIDLEV+1,NLEV,15
          II=II+1
          IC=ICASE(II)
          DO LEV=LV1,MIN(LV1+14,NLEV)
            IC1=IC/4
            ISTEP=IC-4*IC1
            IOC=(ISTEP+1)/2
            SUM=SUM+DBLE(IOC)*ETA(LEV)
            IC=IC1
          END DO
        END DO
        DO J=1,NDWN
          EMU(I,J)=EMU(I,J)+SUM
        END DO
      END DO
C THEN THE LOWER HALF:
      DO I=1,NDWN
        II=NIPWLK*(I-1)
        SUM=0.0D0
        DO LV1=1,MIDLEV,15
        II=II+1
        IC=JCASE(II)
        DO LEV=LV1,MIN(LV1+14,MIDLEV)
          IC1=IC/4
          ISTEP=IC-4*IC1
          IOC=(ISTEP+1)/2
          SUM=SUM+DBLE(IOC)*ETA(LEV)
          IC=IC1
          END DO
        END DO
        DO J=1,NUP
          EMU(J,I)=EMU(J,I)+SUM
        END DO
      END DO

      RETURN
      END
