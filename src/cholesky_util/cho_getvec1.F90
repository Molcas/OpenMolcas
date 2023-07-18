!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,           &
     &                       SCR,LSCR)
!
!     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
!              of symmetry ISYM from file. The vectors are returned
!              in the "current" reduced set. This routine attempts
!              to minimize gather/scatter operations and uses batched
!              reading to (hopefully) improve buffering.
!
!     NOTE: the scratch array SCR(LSCR) is used to read vectors from
!           disk and should not be smaller than NNBSTR(ISYM,1)+1.
!
      use ChoArr, only: iScr
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      REAL*8 CHOVEC(LENVEC,NUMVEC)
      REAL*8 SCR(LSCR)
#include "cholesky.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1)

!     Some initializations.
!     ---------------------

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

!     Get reduced sets of first and last vector.
!     ------------------------------------------

      IRED1 = INFVEC(IVEC1,2,ISYM)
      IRED2 = INFVEC(IVEC2,2,ISYM)

!     Loop through reduced sets to be read.
!     -------------------------------------

      KVEC1 = 1
      JVEC1 = IVEC1
      DO IRED = IRED1,IRED2

!        Count vectors in this reduced set.
!        ----------------------------------

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

!        Skip if this reduced set is empty.
!        ----------------------------------

         IF (JNUM .EQ. 0) GO TO 100

!        Check vector range.
!        -------------------

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

!        Read reduced set index arrays.
!        ------------------------------

         CALL CHO_GETRED(IRED,ILOC,.FALSE.)
         CALL CHO_SETREDIND(ILOC)

!        If reduced sets are identical, simply read the vectors
!        directly into CHOVEC array and go to next reduced set.
!        ------------------------------------------------------

         IF (NNBSTR(ISYM,3) .EQ. NNBSTR(ISYM,2)) THEN
!           IF (CHO_ADRVEC .EQ. 1) THEN
!              IOPT = 2
!              IADR = INFVEC(JVEC1,3,ISYM)
!              LTOT = NNBSTR(ISYM,2)*JNUM
!              CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1),LTOT,IADR)
!              NSYS_CALL = NSYS_CALL + 1
!           ELSE IF (CHO_ADRVEC .EQ. 2) THEN
!              IOPT = 2
!              LTOT = NNBSTR(ISYM,2)
!              DO KK = 1,JNUM
!                 IADR = INFVEC(JVEC1+KK-1,3,ISYM)
!                 CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1+KK-1),
!    &                         LTOT,IADR)
!                 NSYS_CALL = NSYS_CALL + 1
!              END DO
!           ELSE
!              CALL CHO_QUIT('[1] CHO_ADRVEC error in '//SECNAM,102)
!           END IF
!-tbp: replaced above to make use of buffer in cho_vecrd.
            LTOT = NNBSTR(ISYM,2)*JNUM
            JVEC_END = JVEC1 + JNUM - 1
            JNUM_RD = 0
            IREDC = IRED
            MUSED = 0
            CALL CHO_VECRD(CHOVEC(1,KVEC1),LTOT,JVEC1,JVEC_END,ISYM,    &
     &                     JNUM_RD,IREDC,MUSED)
            IF (JNUM_RD .NE. JNUM) THEN
               CALL CHO_QUIT('Logical error [RD1] in '//SECNAM,103)
            END IF
            NSYS_CALL = NSYS_CALL + 1
            GO TO 100
         END IF

!        Set up batch over vectors to be read.
!        -------------------------------------

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

!        Set up mapping between reduced sets.
!        ------------------------------------

         CALL CHO_RS2RS(ISCR,SIZE(ISCR),2,3,IRED,ISYM)

!        Start batch loop.
!        -----------------

         DO IBATCH = 1,NBATCH

            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = JNUM - NVEC*(NBATCH - 1)
            ELSE
               NUMV = NVEC
            END IF
            IBVEC1 = JVEC1 + NVEC*(IBATCH - 1)
            KBVEC1 = KVEC1 + NVEC*(IBATCH - 1)

!           Read vectors.
!           -------------

!           IF (CHO_ADRVEC .EQ. 1) THEN
!              IOPT = 2
!              LTOT = NNBSTR(ISYM,3)*NUMV
!              IADR = INFVEC(IBVEC1,3,ISYM)
!              CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LTOT,IADR)
!              NSYS_CALL = NSYS_CALL + 1
!           ELSE IF (CHO_ADRVEC .EQ. 2) THEN
!              IOPT = 2
!              LTOT = NNBSTR(ISYM,3)
!              KTRG = KSCR
!              DO KK = 1,NUMV
!                 IADR = INFVEC(IBVEC1+KK-1,3,ISYM)
!                 CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KTRG),LTOT,IADR)
!                 KTRG = KTRG + NNBSTR(ISYM,3)
!                 NSYS_CALL = NSYS_CALL + 1
!              END DO
!           ELSE
!              CALL CHO_QUIT('[2] CHO_ADRVEC error in '//SECNAM,102)
!           END IF
!-tbp: replaced above to make use of buffer in cho_vecrd.
            LTOT = NNBSTR(ISYM,3)*NUMV
            JVEC_END = IBVEC1 + NUMV - 1
            JNUM_RD = 0
            IREDC = IRED
            MUSED = 0
            CALL CHO_VECRD(SCR(KSCR),LTOT,IBVEC1,JVEC_END,ISYM,         &
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

!           Copy vectors into result array.
!           -------------------------------

            DO JVEC = 1,NUMV
               KVEC = KBVEC1 + JVEC - 1
               IOFF(1) = IOFF(0) + NNBSTR(ISYM,3)*(JVEC - 1)
               DO IAB = 1,NNBSTR(ISYM,2)
                  KOFF = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                  CHOVEC(IAB,KVEC) = SCR(KOFF)
               END DO
            END DO

         END DO

!        Set next vector to be treated.
!        ------------------------------

  100    KVEC1 = KVEC1 + JNUM
         JVEC1 = JVEC1 + JNUM

      END DO

      END
