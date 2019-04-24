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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      FUNCTION NCASTR_MCLR(IAC,NSTTPI,NTPSTI,ICLSI,NOBATP,NOBTP,IELPTP)
*
* Number of allowed annihilation/creations from a given group
* of stings
*
*     Jeppe Olsen, June 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  NSTTPI(NTPSTI)
      INTEGER  NOBATP(NOBTP)
      INTEGEr IELPTP(NOBTP,*)
*
      LCA = 0
C      write(6,*) ' NCASTR: IELPTP array'
C      write(6,*) ' ===================='
C      write(6,*) ' Number of occupation types ', NTPSTI
C      call iwrtma(IELPTP,NOBTP,NTPSTI,NOBTP,NTPSTI)
C      write(6,*) ' NOBTP NTPSTI ',NOBTP,NTPSTI
C      write(6,*) ' Number of strings per type:'
C      write(6,*) (NSTTPI(I),I=1,NTPSTI)
      DO IOBTP = 1, NOBTP
        DO ISTTP = 1, NTPSTI
*. Type of resulting string
          CALL NEWTYP_MCLR(ICLSI,ISTTP,[IAC],[IOBTP],1,ICLSO,ITPO)
C?        WRITE(6,*) ' IOBTP ISTTP => ITPO, ICLSO '
C?        WRITE(6,*)   IOBTP,ISTTP,ITPO,ICLSO
C          write(6,*) ' IELPTP = ',IELPTP(IOBTP,ISTTP)
C          write(6,*) ' NOBATP = ',NOBATP(IOBTP)
          IF(IAC.EQ.1) THEN
            NENTRY = IELPTP(IOBTP,ISTTP)
          ELSE
            NENTRY = NOBATP(IOBTP)-IELPTP(IOBTP,ISTTP)
          END IF
C          write(6,*) ' NENTRY = ',NENTRY
          IF(ITPO.GT.0) THEN
            LCA = LCA + NENTRY*NSTTPI(ISTTP)
          END IF
*
        END DO
      END DO
*
*     WRITE(6,*) ' Number of generated strings ', LCA
      NCASTR_MCLR = LCA
*
      RETURN
      END
