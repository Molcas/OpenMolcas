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
      SUBROUTINE MKMAW(IDOWN,IDAW,IUP,IRAW,IMAW)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "pt2_guga.fh"
      DIMENSION IDOWN(NVERT,0:3),IDAW(NVERT,0:4)
      DIMENSION IUP(NVERT,0:3),IRAW(NVERT,0:4)
      DIMENSION IMAW(NVERT,0:3)

C COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO IMAW:
      DO 210 IV=MIDV1,NVERT
        DO 211 IC=0,3
          IMAW(IV,IC)=IDAW(IV,IC)
  211   CONTINUE
  210 CONTINUE
C COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO IMAW. HOWEVER,
C    NOTE THAT THE IMAW TABLE IS ACCESSED BY THE UPPER VERTEX.
      DO 230 IU=1,MIDV1-1
        DO 220 IC=0,3
          ID=IDOWN(IU,IC)
          IMAW(IU,IC)=0
          IF(ID.NE.0) IMAW(IU,IC)=IRAW(ID,IC)
  220   CONTINUE
  230 CONTINUE
C FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
      ISUM=1
      DO 250 IV=MIDV1,MIDV2
        DO 240 IC=0,3
          IU=IUP(IV,IC)
          IF(IU.EQ.0) GOTO 240
          IMAW(IU,IC)=ISUM+IMAW(IU,IC)
 240    CONTINUE
        ISUM=ISUM+IRAW(IV,4)
 250  CONTINUE
      DO 270 IV=MIDV1,MIDV2
        DO 260 IC=0,3
          IF(IDOWN(IV,IC).EQ.0) GOTO 260
          IMAW(IV,IC)=ISUM+IMAW(IV,IC)
 260    CONTINUE
        ISUM=ISUM+IDAW(IV,4)
 270  CONTINUE
#ifdef _DEBUGPRINT_
 1010 FORMAT(1X,I4,5X,5(1X,I6))
        WRITE(6,*)
        WRITE(6,*)' THE MODIFIED ARC WEIGHT TABLE IN MKMAW:'
        DO 280 IV=1,NVERT
          WRITE(6,1010) IV,(IMAW(IV,IC),IC=0,3)
  280   CONTINUE
        WRITE(6,*)
#endif
      END SUBROUTINE MKMAW
