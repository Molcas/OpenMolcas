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
      SUBROUTINE CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       SCR,LSCR)
!
!     Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
!              of symmetry ISYM from file. The vectors are returned
!              in the "current" reduced set. This routine attempts
!              to minimize gather/scatter operations along the lines
!              of cho_getvec1. However, in this version, buffering
!              is (hopefully) further improved by reading vectors
!              across reduced sets.
!
!     NOTE: the scratch array SCR(LSCR) is used to read vectors from
!           disk and should not be smaller than NNBSTR(ISYM,1)+1,
!           preferably more.
!
      use ChoArr, only: iScr
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      REAL*8 SCR(LSCR)
#include "cholesky.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC2')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOFF(0:1)

!     Some initializations.
!     ---------------------

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

!     Start buffer batch loop.
!     ------------------------

      KVEC1 = 1
      JVEC1 = IVEC1
      IREDC = -1
      IMAPC = -1
      DO WHILE (JVEC1 .LE. IVEC2)

!        Read as many vectors as fit into scratch space.
!        -----------------------------------------------

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

!        Quit if no vectors were read.
!        -----------------------------

         IF (NVRD .LT. 1) THEN
            CALL CHO_QUIT('Insufficient scratch space for read in '
     &                    //SECNAM,101)
         END IF

!        Loop over reduced sets in scratch space.
!        ----------------------------------------

         JVEC2   = JVEC1 + NVRD - 1
         JRED2   = INFVEC(JVEC2,2,ISYM)
         LVEC1   = JVEC1
         IOFF(1) = KSCR - 1
         DO JRED = JRED1,JRED2

!           Count vectors read from this reduced set.
!           -----------------------------------------

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

!              Read index arrays for this reduced set (if needed).
!              ---------------------------------------------------

               IF (JRED .NE. IREDC) THEN
                  CALL CHO_GETRED(JRED,ILOC,.FALSE.)
                  CALL CHO_SETREDIND(ILOC)
                  IREDC = JRED
               END IF

!              Set up rs-to-rs map (if needed).
!              --------------------------------

               IF (JRED .NE. IMAPC) THEN
                  CALL CHO_RS2RS(ISCR,SIZE(ISCR),2,3,JRED,ISYM)
                  IMAPC = JRED
               END IF

!              Copy vectors to result array.
!              -----------------------------

               DO LVEC = 1,LNUM
                  KVEC = KVEC1 + LVEC - 1
                  DO IAB = 1,NNBSTR(ISYM,2)
                     KOFF = IOFF(MIN(ISCR(IAB),1)) + ISCR(IAB)
                     CHOVEC(IAB,KVEC) = SCR(KOFF)
                  END DO
                  IOFF(1) = IOFF(1) + NNBSTR(ISYM,3)
               END DO

!              Update local vector counters.
!              -----------------------------

               KVEC1 = KVEC1 + LNUM
               LVEC1 = LVEC1 + LNUM

            END IF

         END DO

!        Update global vector counter.
!        -----------------------------

         JVEC1 = JVEC1 + NVRD

      END DO

      END
