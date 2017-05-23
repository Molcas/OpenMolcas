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
      SUBROUTINE CHO_MCA_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,
     &                          LUPRI)
C
C     Purpose: get next keyword and convert it to internal IDKEY.
C
#include "implicit.fh"
      CHARACTER*(*) OPTION(NOPTION)  ! <-- character*(loption)

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_GETKEY')

      PARAMETER (LBLINE = 180, LKWORD = 180)
      CHARACTER*(LBLINE) BLINE
      CHARACTER*(LKWORD) KWORD, KEY

      CHARACTER*(LKWORD) GET_LN
      EXTERNAL GET_LN

      CHARACTER*1 COMMENT
      PARAMETER (COMMENT = '*')

      INTEGER  CHO_TABIND
      EXTERNAL CHO_TABIND

      LOGICAL USE_OBS, USE_ALI

      PARAMETER (NTABLE = 58, LKEY = 4, NEOINP = 1)
      CHARACTER*(LKEY) TABLE(NTABLE)
      CHARACTER*(LKEY) EOINP(NEOINP)

      PARAMETER (NOBSOL = 1, NALIAS = 12)
      CHARACTER*(LKEY) OBSOL(NOBSOL)
      CHARACTER*(LKEY) ALIAS(NALIAS,2)

      DATA TABLE /'THRC','PRIN','BUFF','THRD','DMP1',
     &            'DMP2','SPAN','MINQ','MXSH','SCRE',
     &            'NOSC','QUAL','THRN','WARN','TOON',
     &            'CHEC','CHKA','RSTD','RSTC','RSTM',
     &            'MAXQ','CHOM','REDM','CHKS','CHKM',
     &            'ABSO','NOAB','TRCN','IOVE','REOR',
     &            'HALT','FRAC','QFRA','MXSU','ADDR',
     &            'IFCS','ONES','TWOS','NAIV','VBUF',
     &            'DIAC','TSTS','SSCR','NOSS','SSTH',
     &            'SSNO','1-CE','NO2-','NOPR','PRES',
     &            'PRET','PARA','SIMP','SIMR','FAKE',
     &            'TRUE','BLOC','IDLE'/
      DATA EOINP /'ENDC'/
      DATA OBSOL /'XXXX'/
      DATA ALIAS /'PREC','THSI','THSU','STOP','MEMQ',
     &            'IOMO','DYNA','1CEN','NO2C','THRP',
     &            '1CCD','1C-C',
     &            'THRC','DMP1','DMP2','HALT','QFRA',
     &            'ADDR','VBUF','1-CE','NO2-','PRET',
     &            '1-CE','1-CE'/

C Set flags for using obsolete/alias keywords:
      PARAMETER (USE_OBS = .FALSE., USE_ALI = .TRUE.)

C     Check that we are in sync with caller.
C     --------------------------------------

      IF (NOPTION .NE. NTABLE) THEN
         WRITE(LUPRI,*) SECNAM,': NOPTION = ',NOPTION,
     &                          ' NTABLE = ',NTABLE
         IDKEY = -5
         GO TO 2000
      END IF

C     Other internal checks.
C     ----------------------

      IF (LKEY .GT. LKWORD) THEN
         WRITE(LUPRI,*) SECNAM,': LKEY = ',LKEY,
     &                          ' LKWORD = ',LKWORD
         IDKEY = -5
         GO TO 2000
      END IF

C     Set blank line.
C     ---------------

      DO I = 1,LBLINE
         BLINE(I:I) = ' '
      END DO

C     Read keyword (the check for comment/blank line should be
C     obsolete).
C     --------------------------------------------------------

    1 KEY = GET_LN(LUNIT)
      KWORD = KEY
      CALL UPCASE(KWORD)
      CALL LEFTAD(KWORD)
      DO WHILE (KWORD(1:1).EQ.COMMENT .OR. KWORD.EQ.BLINE)
         KEY = GET_LN(LUNIT)
         KWORD = KEY
         CALL UPCASE(KWORD)
         CALL LEFTAD(KWORD)
      END DO

      LAST = ICLAST(KWORD,LEN(KWORD))
      DO I = LAST+1,LKEY
         KWORD(I:I) = ' '
      END DO

C     Check for obsolete keyword.
C     ---------------------------

      IF (USE_OBS) THEN
         IOBSOL = CHO_TABIND(OBSOL,LKEY,NOBSOL,' ',0,0,KWORD(1:LKEY))
         IF (IOBSOL.GT.0 .AND. IOBSOL.LE.NOBSOL) THEN
            WRITE(LUPRI,*) '*** NOTICE: Cholesky keyword "',
     &      KWORD(1:LKEY),'" is obsolete and will be disregarded.'
            GO TO 1
         END IF
      END IF

C     Check for alias.
C     ----------------

      IALIAS = 0
      IF (USE_ALI) THEN
         IALIAS = CHO_TABIND(ALIAS(1,1),LKEY,NALIAS,' ',0,0,
     &                       KWORD(1:LKEY))
         IF (IALIAS.GT.0 .AND. IALIAS.LE.NALIAS) THEN
            KWORD(1:LKEY) = ALIAS(IALIAS,2)
         ELSE
            IALIAS = 0
         END IF
      END IF

C     Table lookup.
C     -------------

      IDKEY = CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LKEY,NEOINP,
     &                   KWORD(1:LKEY))
      IF (IDKEY .EQ. -1) THEN
         WRITE(LUPRI,*) SECNAM,': keyword not recognized:'
         IF (IALIAS .GT. 0) THEN
            WRITE(LUPRI,*) 'Internal  key: ',KWORD(1:LAST),
     &                     ' (significant part: ',KWORD(1:LKEY),')'
            WRITE(LUPRI,*) 'Aliasing used: ',ALIAS(IALIAS,1),
     &                     ' <-> ',ALIAS(IALIAS,2)
         ELSE
            WRITE(LUPRI,*) 'Internal  key: ',KWORD(1:LAST),
     &                     ' (significant part: ',KWORD(1:LKEY),')'
         END IF
         WRITE(LUPRI,*)
         IF (LOPTION .GT. 0) THEN
            WRITE(LUPRI,*) 'Available keywords and short explanations:'
            DO IOPTION = 1,NOPTION
               IF (TABLE(IOPTION) .NE. 'XXXX') THEN
                  WRITE(LUPRI,*) TABLE(IOPTION),': ',
     &                           OPTION(IOPTION)
               END IF
            END DO
         ELSE
            WRITE(LUPRI,*) 'Available keywords:'
            DO IOPTION = 1,NOPTION
               IF (TABLE(IOPTION) .NE. 'XXXX') THEN
                  WRITE(LUPRI,*) TABLE(IOPTION)
               END IF
            END DO
         END IF
         WRITE(LUPRI,*)
      END IF

C     Normal exit point.
C     ------------------

 2000 RETURN

      END
