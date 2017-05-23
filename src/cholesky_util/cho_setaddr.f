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
      SUBROUTINE CHO_SETADDR(INFRED,INFVEC,MRED,MVEC,M2,MSYM)
C
C     Purpose: set first disk addresses for reduced set info and
C              vectors.
C
#include "implicit.fh"
      INTEGER INFRED(MRED), INFVEC(MVEC,M2,MSYM)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_SETADDR')

C     Set addresses.
C     --------------

      IF (XNPASS .EQ. 0) THEN
         INFRED(1) = 0
         DO ISYM = 1,NSYM
            INFVEC(1,3,ISYM) = 0
            INFVEC(1,4,ISYM) = 0
         END DO
      ELSE IF (XNPASS .GT. 0) THEN
         IRED  = 3
         IPASS = XNPASS
         KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(IRED - 1)
         KOFF2 = ip_INDRED   + MMBSTRT*(IRED - 1)
         CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                   IWORK(KOFF2),IWORK(ip_INDRSH),IWORK(ip_iSP2F),
     &                   MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,.FALSE.)
         CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),IWORK(ip_NNBSTRSH),
     &                      NSYM,NNSHL,IRED)
         IF (IPASS .EQ. 1) THEN
            INFRED(IPASS+1) = INFRED(IPASS)
     &                      + NSYM*NNSHL + 2*NNBSTRT(IRED) + NNSHL
         ELSE
            INFRED(IPASS+1) = INFRED(IPASS)
     &                      + NSYM*NNSHL + NNBSTRT(IRED)
         END IF
         DO ISYM = 1,NSYM
            IF (NUMCHO(ISYM) .EQ. 0) THEN
               INFVEC(1,3,ISYM) = 0
               INFVEC(1,4,ISYM) = 0
            ELSE IF (NUMCHO(ISYM) .GT. 0) THEN
               IF (CHO_ADRVEC .EQ. 1) THEN
                  JPASS = INFVEC(NUMCHO(ISYM),2,ISYM)
                  IF (JPASS .EQ. IPASS) THEN
                     INFVEC(NUMCHO(ISYM)+1,3,ISYM) =
     &               INFVEC(NUMCHO(ISYM),3,ISYM) + NNBSTR(ISYM,IRED)
                     INFVEC(NUMCHO(ISYM)+1,4,ISYM) =
     &               INFVEC(NUMCHO(ISYM),4,ISYM) + NNBSTR(ISYM,IRED)
                  ELSE IF (JPASS.LE.XNPASS .AND. JPASS.GT.0) THEN
                     IPASS = JPASS
                     KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(IRED - 1)
                     KOFF2 = ip_INDRED   + MMBSTRT*(IRED - 1)
                     CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                               IWORK(KOFF2),IWORK(ip_INDRSH),
     &                               IWORK(ip_iSP2F),
     &                               MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,
     &                               .FALSE.)
                     CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),
     &                                  IWORK(ip_NNBSTRSH),
     &                                  NSYM,NNSHL,IRED)
                     INFVEC(NUMCHO(ISYM)+1,3,ISYM) =
     &               INFVEC(NUMCHO(ISYM),3,ISYM) + NNBSTR(ISYM,IRED)
                     INFVEC(NUMCHO(ISYM)+1,4,ISYM) =
     &               INFVEC(NUMCHO(ISYM),4,ISYM) + NNBSTR(ISYM,IRED)
                  ELSE
                     CALL CHO_QUIT('[1] JPASS error in '//SECNAM,104)
                  END IF
               ELSE IF (CHO_ADRVEC .EQ. 2) THEN
                  JPASS = INFVEC(NUMCHO(ISYM),2,ISYM)
                  IF (JPASS .EQ. IPASS) THEN
                     LSA = NNBSTR(ISYM,IRED)
                     CALL CHO_MEM('SetAddr','ALLO','REAL',KSA,LSA)
                     IOPT = 2
                     IADR = INFVEC(NUMCHO(ISYM),3,ISYM)
                     CALL DDAFILE(LUCHO(ISYM),IOPT,WORK(KSA),LSA,IADR)
                     INFVEC(NUMCHO(ISYM)+1,3,ISYM) = IADR
                     INFVEC(NUMCHO(ISYM)+1,4,ISYM) =
     &                 INFVEC(NUMCHO(ISYM),4,ISYM) + NNBSTR(ISYM,IRED)
                     CALL CHO_MEM('SetAddr','FREE','REAL',KSA,LSA)
                  ELSE IF (JPASS.LE.XNPASS .AND. JPASS.GT.0) THEN
                     IPASS = JPASS
                     KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(IRED - 1)
                     KOFF2 = ip_INDRED   + MMBSTRT*(IRED - 1)
                     CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                               IWORK(KOFF2),IWORK(ip_INDRSH),
     &                               IWORK(ip_iSP2F),
     &                               MAXRED,NSYM,NNSHL,MMBSTRT,IPASS,
     &                               .FALSE.)
                     CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),
     &                                  IWORK(ip_NNBSTRSH),
     &                                  NSYM,NNSHL,IRED)
                     LSA = NNBSTR(ISYM,IRED)
                     CALL CHO_MEM('SetAddr','ALLO','REAL',KSA,LSA)
                     IOPT = 2
                     IADR = INFVEC(NUMCHO(ISYM),3,ISYM)
                     CALL DDAFILE(LUCHO(ISYM),IOPT,WORK(KSA),LSA,IADR)
                     INFVEC(NUMCHO(ISYM)+1,3,ISYM) = IADR
                     INFVEC(NUMCHO(ISYM)+1,4,ISYM) =
     &                 INFVEC(NUMCHO(ISYM),4,ISYM) + NNBSTR(ISYM,IRED)
                     CALL CHO_MEM('SetAddr','FREE','REAL',KSA,LSA)
                  ELSE
                     CALL CHO_QUIT('[2] JPASS error in '//SECNAM,104)
                  END IF
               ELSE
                  CALL CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
               END IF
            ELSE
               CALL CHO_QUIT('NUMCHO error in '//SECNAM,104)
            END IF
         END DO
      ELSE
         CALL CHO_QUIT('XNPASS error in '//SECNAM,104)
      END IF

      END
