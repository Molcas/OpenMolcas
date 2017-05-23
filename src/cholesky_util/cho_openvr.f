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
      SUBROUTINE CHO_OPENVR(IOPT,ID)
C
C     Purpose: open (IOPT=1) or close (IOPT=2) files for vector
C              and reduced set storage as well as restart files.
C              ID=1: open local files (for parallel run)
C              ID=2: open global files (as in serial run)
C
#include "implicit.fh"
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_OPENVR')

      CHARACTER*5 FNRED
      CHARACTER*6 FNVEC(8), FRST, FMAP

      IF (IOPT .EQ. 1) THEN
         FMAP = 'CHOMAP'
         IF (ID .EQ. 1) THEN
            FNRED = 'CHRDL'
            DO ISYM = 1,NSYM
               WRITE(FNVEC(ISYM),'(A5,I1)') 'CHVCL',ISYM
            END DO
            FRST = 'CHRSTL'
         ELSE
            FNRED = 'CHRED'
            DO ISYM = 1,NSYM
               WRITE(FNVEC(ISYM),'(A5,I1)') 'CHVEC',ISYM
            END DO
            FRST = 'CHORST'
         END IF
         LURED = 7
         CALL DANAME_MF_WA(LURED,FNRED)
         IF (CHO_ADRVEC .EQ. 1) THEN
            DO ISYM = 1,NSYM
               LUCHO(ISYM) = 7
               CALL DANAME_MF_WA(LUCHO(ISYM),FNVEC(ISYM))
            END DO
         ELSE IF (CHO_ADRVEC .EQ. 2) THEN
            DO ISYM = 1,NSYM
               LUCHO(ISYM) = 7
               CALL DANAME_MF(LUCHO(ISYM),FNVEC(ISYM))
            END DO
         ELSE
            CALL CHO_QUIT('CHO_ADRVEC out of bounds in '//SECNAM,102)
         END IF
         LURST = 7
         CALL DANAME_MF_WA(LURST,FRST)
         LUMAP = 7
         CALL DANAME(LUMAP,FMAP)
      ELSE IF (IOPT .EQ. 2) THEN
         IF (LURED .GT. 0) THEN
            CALL DACLOS(LURED)
            LURED = 0
         END IF
         DO ISYM = 1,NSYM
            IF (LUCHO(ISYM) .GT. 0) THEN
               CALL DACLOS(LUCHO(ISYM))
               LUCHO(ISYM) = 0
            END IF
         END DO
         IF (LURST .GT. 0) THEN
            CALL DACLOS(LURST)
            LURST = 0
         END IF
         IF (LUMAP .GT. 0) THEN
            CALL DACLOS(LUMAP)
            LUMAP = 0
         END IF
      ELSE
         WRITE(LUPRI,*) SECNAM,': IOPT out of bounds: ',IOPT
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

      END
