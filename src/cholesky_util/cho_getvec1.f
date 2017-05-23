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
      SUBROUTINE CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       SCR,LSCR)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set. This routine attempts
C              to minimize gather/scatter operations and uses batched
C              reading to (hopefully) improve buffering.
C
C     NOTE: the scratch array SCR(LSCR) is used to read vectors from
C           disk and should not be smaller than NNBSTR(ISYM,1)+1.
C
#include "implicit.fh"
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      DIMENSION SCR(LSCR)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1)

      PARAMETER (N2 = INFVEC_N2)

      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      ISCR(I)=IWORK(ip_ISCR-1+I)

C     Some initializations.
C     ---------------------

      ILOC  = 3
      IVEC2 = IVEC1 + NUMVEC - 1

      KJUNK = 1
      KSCR  = KJUNK + 1
      LEFT  = LSCR  - KSCR + 1
      IF (LEFT .LE. 0) THEN
         CALL CHO_QUIT('Insufficient scratch space in '//SECNAM,101)
      END IF

      SCR(KJUNK) = 0.0D0
      IOFF(0)    = KJUNK

C     Get reduced sets of first and last vector.
C     ------------------------------------------

      IRED1 = INFVEC(IVEC1,2,ISYM)
      IRED2 = INFVEC(IVEC2,2,ISYM)

C     Loop through reduced sets to be read.
C     -------------------------------------

      KVEC1 = 1
      JVEC1 = IVEC1
      DO IRED = IRED1,IRED2

C        Count vectors in this reduced set.
C        ----------------------------------

         JNUM = 0
         JVEC = JVEC1 - 1
         DO WHILE (JVEC.LT.IVEC2)
            JVEC = JVEC + 1
            JRED = INFVEC(JVEC,2,ISYM)
            IF (JRED .EQ. IRED) THEN
               JNUM = JNUM + 1 ! increase counter
            ELSE
               JVEC = IVEC2 ! break while loop
            END IF
         END DO

C        Skip if this reduced set is empty.
C        ----------------------------------

         IF (JNUM .EQ. 0) GO TO 100

C        Check vector range.
C        -------------------

         IF (LOCDBG) THEN
            JVEC2 = JVEC1 + JNUM - 1
            IF (JVEC2 .GT. IVEC2) THEN
               WRITE(LUPRI,*) SECNAM,': IRED  = ',IRED
               WRITE(LUPRI,*) SECNAM,': JNUM  = ',JNUM
               WRITE(LUPRI,*) SECNAM,': JVEC1 = ',JVEC1
               WRITE(LUPRI,*) SECNAM,': JVEC2 = ',JVEC2
               CALL CHO_QUIT('Vector index error in '//SECNAM,103)
            END IF
         END IF

C        Read reduced set index arrays.
C        ------------------------------

         KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(ILOC - 1)
         KOFF2 = ip_INDRED   + MMBSTRT*(ILOC - 1)
         CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                   IWORK(KOFF2),IWORK(ip_INDRSH),IWORK(ip_iSP2F),
     &                   MAXRED,NSYM,NNSHL,MMBSTRT,IRED,
     &                   .FALSE.)
         CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),
     &                      IWORK(ip_NNBSTRSH),NSYM,NNSHL,3)

C        If reduced sets are identical, simply read the vectors
C        directly into CHOVEC array and go to next reduced set.
C        ------------------------------------------------------

         IF (NNBSTR(ISYM,3) .EQ. NNBSTR(ISYM,2)) THEN
C           IF (CHO_ADRVEC .EQ. 1) THEN
C              IOPT = 2
C              IADR = INFVEC(JVEC1,3,ISYM)
C              LTOT = NNBSTR(ISYM,2)*JNUM
C              CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1),LTOT,IADR)
C              NSYS_CALL = NSYS_CALL + 1
C           ELSE IF (CHO_ADRVEC .EQ. 2) THEN
C              IOPT = 2
C              LTOT = NNBSTR(ISYM,2)
C              DO KK = 1,JNUM
C                 IADR = INFVEC(JVEC1+KK-1,3,ISYM)
C                 CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1+KK-1),
C    &                         LTOT,IADR)
C                 NSYS_CALL = NSYS_CALL + 1
C              END DO
C           ELSE
C              CALL CHO_QUIT('[1] CHO_ADRVEC error in '//SECNAM,102)
C           END IF
C-tbp: replaced above to make use of buffer in cho_vecrd.
            LTOT = NNBSTR(ISYM,2)*JNUM
            JVEC_END = JVEC1 + JNUM - 1
            JNUM_RD = 0
            IREDC = IRED
            MUSED = 0
            CALL CHO_VECRD(CHOVEC(1,KVEC1),LTOT,JVEC1,JVEC_END,ISYM,
     &                     JNUM_RD,IREDC,MUSED)
            IF (JNUM_RD .NE. JNUM) THEN
               CALL CHO_QUIT('Logical error [RD1] in '//SECNAM,103)
            END IF
            NSYS_CALL = NSYS_CALL + 1
            GO TO 100
         END IF

C        Set up batch over vectors to be read.
C        -------------------------------------

         MINL = NNBSTR(ISYM,3)
         IF (MINL .LT. 1) THEN
            NVEC = 0
         ELSE
            NVEC = MIN(LEFT/MINL,JNUM)
         END IF
         IF (NVEC .LT. 1) THEN
            WRITE(LUPRI,*) SECNAM,': insufficient scratch space:'
            WRITE(LUPRI,*) 'LEFT = ',LEFT
            WRITE(LUPRI,*) 'JNUM = ',JNUM
            WRITE(LUPRI,*) 'MINL = ',MINL
            WRITE(LUPRI,*) 'NVEC = ',NVEC
            WRITE(LUPRI,*) 'Input:'
            WRITE(LUPRI,*) 'IVEC1  = ',IVEC1
            WRITE(LUPRI,*) 'NUMVEC = ',NUMVEC
            WRITE(LUPRI,*) 'LENVEC = ',LENVEC
            WRITE(LUPRI,*) 'ISYM   = ',ISYM
            CALL CHO_QUIT('Insufficient scratch space in '//SECNAM,104)
            NBATCH = 0 ! to avoid compiler warnings
         ELSE
            NBATCH = (JNUM - 1)/NVEC + 1
         END IF

C        Set up mapping between reduced sets.
C        ------------------------------------

         CALL CHO_RS2RS(IWORK(ip_ISCR),l_ISCR,2,3,IRED,ISYM)

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = JNUM - NVEC*(NBATCH - 1)
            ELSE
               NUMV = NVEC
            END IF
            IBVEC1 = JVEC1 + NVEC*(IBATCH - 1)
            KBVEC1 = KVEC1 + NVEC*(IBATCH - 1)

C           Read vectors.
C           -------------

C           IF (CHO_ADRVEC .EQ. 1) THEN
C              IOPT = 2
C              LTOT = NNBSTR(ISYM,3)*NUMV
C              IADR = INFVEC(IBVEC1,3,ISYM)
C              CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LTOT,IADR)
C              NSYS_CALL = NSYS_CALL + 1
C           ELSE IF (CHO_ADRVEC .EQ. 2) THEN
C              IOPT = 2
C              LTOT = NNBSTR(ISYM,3)
C              KTRG = KSCR
C              DO KK = 1,NUMV
C                 IADR = INFVEC(IBVEC1+KK-1,3,ISYM)
C                 CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KTRG),LTOT,IADR)
C                 KTRG = KTRG + NNBSTR(ISYM,3)
C                 NSYS_CALL = NSYS_CALL + 1
C              END DO
C           ELSE
C              CALL CHO_QUIT('[2] CHO_ADRVEC error in '//SECNAM,102)
C           END IF
C-tbp: replaced above to make use of buffer in cho_vecrd.
            LTOT = NNBSTR(ISYM,3)*NUMV
            JVEC_END = IBVEC1 + NUMV - 1
            JNUM_RD = 0
            IREDC = IRED
            MUSED = 0
            CALL CHO_VECRD(SCR(KSCR),LTOT,IBVEC1,JVEC_END,ISYM,
     &                     JNUM_RD,IREDC,MUSED)
            IF (JNUM_RD .NE. NUMV) THEN
               CALL CHO_QUIT('Logical error [RD2] in '//SECNAM,103)
            END IF
            IF (CHO_ADRVEC .EQ. 1) THEN
               NSYS_CALL = NSYS_CALL + 1
            ELSE IF (CHO_ADRVEC .EQ. 2) THEN
               NSYS_CALL = NSYS_CALL + NUMV
            ELSE
               CALL CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
            END IF

C           Copy vectors into result array.
C           -------------------------------

            DO JVEC = 1,NUMV
               KVEC = KBVEC1 + JVEC - 1
               IOFF(1) = IOFF(0) + NNBSTR(ISYM,3)*(JVEC - 1)
               DO IAB = 1,NNBSTR(ISYM,2)
                  KOFF = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                  CHOVEC(IAB,KVEC) = SCR(KOFF)
               END DO
            END DO

         END DO

C        Set next vector to be treated.
C        ------------------------------

  100    KVEC1 = KVEC1 + JNUM
         JVEC1 = JVEC1 + JNUM

      END DO

      END
