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
      SUBROUTINE CHO_RS2RS(IMAP,LMAP,IRS2,IRS3,IRED3,ISYM)
C
C     Purpose: set up mapping between reduced sets stored at IRS2 and
C              IRS3 (IRED3 is the reduced set id of IRS3).
C
C     WARNING: for IRED3 = 1, INDRED is reset here!!!!
C
      use ChoSwp, only: nnBstRSh, iiBstRSh
#include "implicit.fh"
      INTEGER IMAP(LMAP)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_RS2RS')

      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)

C     Check input.
C     ------------

      IF (IRS2.LT.1 .OR. IRS2.GT.3 .OR.
     &    IRS3.LT.1 .OR. IRS3.GT.3) THEN
         CALL CHO_QUIT('Index error in '//SECNAM,104)
      ELSE IF (LMAP .LT. NNBSTR(ISYM,IRS2)) THEN
         CALL CHO_QUIT('Dimension error in '//SECNAM,104)
      END IF

C     For IRED3 = 1, INDRED array addresses into shell pair. We hence
C     need to reset it (as warned about above).
C     ---------------------------------------------------------------

      IF (IRED3 .EQ. 1) THEN
         K0 = ip_INDRED - 1 + MMBSTRT*(IRS3 - 1)
         I1 = IIBSTR(ISYM,IRS3) + 1
         I2 = I1 + NNBSTR(ISYM,IRS3) - 1
         DO I = I1,I2
            IWORK(K0+I) = I
         END DO
      END IF

C     Set up mapping array.
C     ---------------------

      CALL CHO_IZERO(IMAP,NNBSTR(ISYM,IRS2))
      DO ISHLAB = 1,NNSHL
         N2 = NNBSTRSH(ISYM,ISHLAB,IRS2)
         N3 = NNBSTRSH(ISYM,ISHLAB,IRS3)
         IF (N2.GT.0 .AND. N3.GT.0) THEN
            IF (N2 .LT. N3) THEN
               IAB1 = IIBSTRSH(ISYM,ISHLAB,IRS2) + 1
               IAB2 = IAB1 + N2 - 1
               LAST = 0
               DO IAB = IAB1,IAB2
                  JAB = INDRED(IIBSTR(ISYM,IRS2)+IAB,IRS2)
                  KAB = LAST
                  DO WHILE (KAB .LT. NNBSTRSH(ISYM,ISHLAB,IRS3))
                     KAB = KAB + 1
                     LAB = IIBSTRSH(ISYM,ISHLAB,IRS3) + KAB
                     MAB = INDRED(IIBSTR(ISYM,IRS3)+LAB,IRS3)
                     IF (MAB .EQ. JAB) THEN
                        IMAP(IAB) = LAB
                        LAST = KAB
                        KAB  = NNBSTRSH(ISYM,ISHLAB,IRS3)
                     END IF
                  END DO
               END DO
            ELSE
               LAB1 = IIBSTRSH(ISYM,ISHLAB,IRS3) + 1
               LAB2 = LAB1 + N3 - 1
               LAST = 0
               DO LAB = LAB1,LAB2
                  MAB = INDRED(IIBSTR(ISYM,IRS3)+LAB,IRS3)
                  KAB = LAST
                  DO WHILE (KAB .LT. NNBSTRSH(ISYM,ISHLAB,IRS2))
                     KAB = KAB + 1
                     IAB = IIBSTRSH(ISYM,ISHLAB,IRS2) + KAB
                     JAB = INDRED(IIBSTR(ISYM,IRS2)+IAB,IRS2)
                     IF (JAB .EQ. MAB) THEN
                        IMAP(IAB) = LAB
                        LAST = KAB
                        KAB  = NNBSTRSH(ISYM,ISHLAB,IRS2)
                     END IF
                  END DO
               END DO
            END IF
         END IF
      END DO

      END
