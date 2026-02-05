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
      use definitions, only: iwp, wp
      use constants, only: Zero
      use gugx, only: SGS, CIS
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)
      integer(kind=iwp), intent(in):: ICASE(*),JCASE(*)
      integer(kind=iwp), intent(in):: NUP,NDWN
      real(kind=wp), intent(out):: EMU(NUP,NDWN)

      Integer(kind=iwp) nLev, nIpWlk
      Integer(kind=iwp) I,II,LV1,IC,LEV,IC1,ISTEP,IOC,J

      nLev  = SGS%nLev
      nIpWlk= CIS%nIpWlk

      DO I=1,NUP
        II=NIPWLK*(I-1)
        SUM=Zero
        DO LV1=SGS%MIDLEV+1,NLEV,15
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
        SUM=Zero
        DO LV1=1,SGS%MIDLEV,15
        II=II+1
        IC=JCASE(II)
        DO LEV=LV1,MIN(LV1+14,SGS%MIDLEV)
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

      END SUBROUTINE DIELMV
