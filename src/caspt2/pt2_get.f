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
      SUBROUTINE PT2_GET(NSIZE,LAB,VEC)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"

      INTEGER NSIZE
      CHARACTER(len=*) LAB
      REAL*8 VEC(*)

      CHARACTER(LEN=8) LAB1

      INTEGER I,IAD,NSZ
#ifdef _DEBUGPRINT_
      INTEGER J
#endif

      I=9-LEN(LAB)
      IF(I.GE.1) THEN
        LAB1='        '
        LAB1(I:8)=LAB
      ELSE
        LAB1=LAB(1:8)
      END IF

C FIND DISK ADDRESS:
      DO I=1,64
        IF (CLAB10(I).EQ.LAB1) THEN
          NSZ=MIN(IADR10(I,2),NSIZE)
          IAD=IADR10(I,1)
          CALL DDAFILE(LUDMAT,2,VEC,NSZ,IAD)
#ifdef _DEBUGPRINT_
          WRITE(6,*) LAB1,' SUCCESSFULLY READ FROM LUDMAT.'
          WRITE(6,*)'         SIZE:',NSZ,' *8 BYTES'
          WRITE(6,*)' DISK ADDRESS:',IADR10(I,1)
          WRITE(6,'(10F12.8)') (VEC(J),J=1,MIN(10,NSZ))
#endif
          RETURN
        END IF
      END DO
      WRITE(6,*)' LABEL ',LAB1,' NOT FOUND ON LUDMAT.'
      CALL ABEND()
      END
