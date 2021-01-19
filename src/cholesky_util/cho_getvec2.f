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
      SUBROUTINE CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       SCR,LSCR)
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set. This routine attempts
C              to minimize gather/scatter operations along the lines
C              of cho_getvec1. However, in this version, buffering
C              is (hopefully) further improved by reading vectors
C              across reduced sets.
C
C     NOTE: the scratch array SCR(LSCR) is used to read vectors from
C           disk and should not be smaller than NNBSTR(ISYM,1)+1,
C           preferably more.
C
      use ChoArr, only: iSP2F, iScr
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, InfRed, InfVec,
     &                  IndRed
#include "implicit.fh"
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      DIMENSION SCR(LSCR)
#include "cholesky.fh"
#include "choptr.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC2')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1)

C     Some initializations.
C     ---------------------

      ILOC = 3

      IVEC2 = IVEC1 + NUMVEC - 1

      KJUNK = 1
      KSCR  = KJUNK + 1
      LEFT  = LSCR  - KSCR + 1
      IF (LEFT .LT. 1) THEN
         CALL CHO_QUIT('Insufficient scratch space in '//SECNAM,101)
      END IF

      SCR(KJUNK) = 0.0D0
      IOFF(0)    = KJUNK

C     Start buffer batch loop.
C     ------------------------

      KVEC1 = 1
      JVEC1 = IVEC1
      IREDC = -1
      IMAPC = -1
      DO WHILE (JVEC1 .LE. IVEC2)

C        Read as many vectors as fit into scratch space.
C        -----------------------------------------------

         JRED1 = INFVEC(JVEC1,2,ISYM)
         NVRD  = 0
         MUSED = 0
         CALL CHO_VECRD(SCR(KSCR),LEFT,JVEC1,IVEC2,ISYM,
     &                  NVRD,IREDC,MUSED)
         IF (CHO_ADRVEC .EQ. 1) THEN
            NSYS_CALL = NSYS_CALL + 1
         ELSE IF (CHO_ADRVEC .EQ. 2) THEN
            NSYS_CALL = NSYS_CALL + NVRD
         ELSE
            CALL CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
         END IF

C        Quit if no vectors were read.
C        -----------------------------

         IF (NVRD .LT. 1) THEN
            CALL CHO_QUIT('Insufficient scratch space for read in '
     &                    //SECNAM,101)
         END IF

C        Loop over reduced sets in scratch space.
C        ----------------------------------------

         JVEC2   = JVEC1 + NVRD - 1
         JRED2   = INFVEC(JVEC2,2,ISYM)
         LVEC1   = JVEC1
         IOFF(1) = KSCR - 1
         DO JRED = JRED1,JRED2

C           Count vectors read from this reduced set.
C           -----------------------------------------

            LNUM = 0
            LVEC = LVEC1 - 1
            DO WHILE (LVEC.LT.JVEC2)
               LVEC = LVEC + 1
               LRED = INFVEC(LVEC,2,ISYM)
               IF (LRED .EQ. JRED) THEN
                  LNUM = LNUM + 1 ! increase counter
               ELSE
                  LVEC = JVEC2 ! break loop
               END IF
            END DO

            IF (LNUM .GT. 0) THEN

C              Read index arrays for this reduced set (if needed).
C              ---------------------------------------------------

               IF (JRED .NE. IREDC) THEN
                  CALL CHO_GETRED(INFRED,nnBstRSh(:,:,ILOC),
     &                            IndRed(1,ILOC),INDRSH,iSP2F,
     &                            MAXRED,NSYM,NNSHL,MMBSTRT,JRED,
     &                            .FALSE.)
                  CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,3)
                  IREDC = JRED
               END IF

C              Set up rs-to-rs map (if needed).
C              --------------------------------

               IF (JRED .NE. IMAPC) THEN
                  CALL CHO_RS2RS(ISCR,SIZE(ISCR),2,3,JRED,ISYM)
                  IMAPC = JRED
               END IF

C              Copy vectors to result array.
C              -----------------------------

               DO LVEC = 1,LNUM
                  KVEC = KVEC1 + LVEC - 1
                  DO IAB = 1,NNBSTR(ISYM,2)
                     KOFF = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                     CHOVEC(IAB,KVEC) = SCR(KOFF)
                  END DO
                  IOFF(1) = IOFF(1) + NNBSTR(ISYM,3)
               END DO

C              Update local vector counters.
C              -----------------------------

               KVEC1 = KVEC1 + LNUM
               LVEC1 = LVEC1 + LNUM

            END IF

         END DO

C        Update global vector counter.
C        -----------------------------

         JVEC1 = JVEC1 + NVRD

      END DO

      END
