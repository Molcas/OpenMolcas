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
      SUBROUTINE CHO_PUTRED1(INFRED,NNBSTRSH,INDRED,INDRSH,ISP2F,
     &                       MRED,MSYM,MMSHL,LMMBSTRT,IPASS,ILOC)
C
C     Purpose: write index arrays for current reduced set (reduced set
C              IPASS).
C
#include "implicit.fh"
      INTEGER INFRED(MRED)
      INTEGER NNBSTRSH(MSYM,MMSHL), INDRED(LMMBSTRT), INDRSH(LMMBSTRT)
      INTEGER ISP2F(MMSHL)
#include "cholesky.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_PUTRED1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Test dimensions.
C     ----------------

      IF (ILOC.LT.1 .OR. ILOC.GT.3) THEN
         CALL CHO_QUIT('ILOC error in '//SECNAM,104)
      END IF

      IF (MSYM .NE. NSYM) THEN
         CALL CHO_QUIT('NSYM error in '//SECNAM,104)
      END IF

      IF (MMSHL .NE. NNSHL) THEN
         CALL CHO_QUIT('NNSHL error in '//SECNAM,104)
      END IF

      IF (LMMBSTRT .NE. NNBSTRT(1)) THEN
         CALL CHO_QUIT('NNBSTRT(1) error in '//SECNAM,104)
      END IF

      IF (LMMBSTRT .LT. NNBSTRT(ILOC)) THEN
         CALL CHO_QUIT('NNBSTRT(ILOC) error in '//SECNAM,104)
      END IF

      IF ((IPASS.LT.1) .OR. (IPASS.GT.MAXRED)) THEN
         CALL CHO_QUIT('IPASS error in '//SECNAM,104)
      END IF

C     Get first address.
C     ------------------

      IADR1 = INFRED(IPASS)
      IF (IADR1 .LT. 0) THEN
         WRITE(LUPRI,*) SECNAM,': negative address for reduced set ',
     &                  IPASS,': ',IADR1
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': putting reduced set ',IPASS,
     &                  ' at addr: ',IADR1
      END IF

C     Write index arrays.
C     -------------------

      IOPT = 1
      IADR = IADR1
      LTOT = NSYM*NNSHL
      CALL IDAFILE(LURED,IOPT,NNBSTRSH(1,1),LTOT,IADR)
      IOPT = 1
      IADR = IADR1 + NSYM*NNSHL
      LTOT = NNBSTRT(ILOC)
      CALL IDAFILE(LURED,IOPT,INDRED,LTOT,IADR)
      IF (IPASS .EQ. 1) THEN
         IOPT = 1
         IADR = IADR1 + NSYM*NNSHL + NNBSTRT(1)
         LTOT = NNBSTRT(1)
         CALL IDAFILE(LURED,IOPT,INDRSH,LTOT,IADR)
         IOPT = 1
         IADR = IADR1 + NSYM*NNSHL + 2*NNBSTRT(1)
         LTOT = NNSHL
         CALL IDAFILE(LURED,IOPT,ISP2F,LTOT,IADR)
      END IF

      END
