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
      SUBROUTINE PT2_PUT(NSIZE,LAB,VEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VEC(*)
      CHARACTER(8) LAB1
      CHARACTER(len=*) LAB
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"

      I=9-LEN(LAB)
      IF(I.GE.1) THEN
        LAB1='        '
        LAB1(I:8)=LAB
      ELSE
        LAB1=LAB(1:8)
      END IF

C FIND DISK ADDRESS:
      DO I=1,64
        IF(CLAB10(I).EQ.'   EMPTY') THEN
          CLAB10(I)=LAB1
          IAD=IADR10(I,1)
          IADR10(I,2)=NSIZE
          CALL DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
          IF(I.LT.64) IADR10(I+1,1)=IAD
          GOTO 20
        ELSE IF (CLAB10(I).EQ.LAB1) THEN
          IF(NSIZE.GT.IADR10(I,2)) GOTO 98
          IAD=IADR10(I,1)
          IADR10(I,2)=NSIZE
          CALL DDAFILE(LUDMAT,1,VEC,NSIZE,IAD)
          GOTO 20
        END IF
      END DO
      WRITE(6,*)' NO MORE AVAILABLE FIELDS ON FILE DENS.'
      GOTO 99
  20  CONTINUE
      RETURN
  98  WRITE(6,*)' ATTEMPT TO INCREASE SIZE OF A FIELD.'
  99  WRITE(6,*)' SUBROUTINE PUT FAILS.'
      CALL ERRTRA
      CALL ABEND()
      END
