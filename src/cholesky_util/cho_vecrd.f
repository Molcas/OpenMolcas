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
      SUBROUTINE CHO_VECRD(SCR,LSCR,JVEC1,IVEC2,ISYM,
     &                     JNUM,IREDC,MUSED)
!
!     Purpose: read as many vectors as fit into SCR array starting
!              at vector JVEC1 and reading at most until vector IVEC2.
!              On exit, JNUM is the number of vectors read.
!              On entry as well as exit, IREDC identifies the reduced
!              set stored in core (at position "3"; use -1 if none
!              or unkown). Vectors are taken from buffer, if possible.
!
!     NOTE: if no vectors can be read, JNUM=0 and MUSED=0 are returned,
!           but execution is NOT stopped here!!!
!
      Implicit Real*8 (a-h,o-z)
      DIMENSION SCR(LSCR)
#include "cholesky.fh"

      LOGICAL DOREAD

!     Initialize.
!     -----------

      JNUM = 0
      MUSED = 0
      IF (LSCR .LT. 1) RETURN

!     Copy vectors from buffer (if possible).
!     Only relevant for "external" runs (else the vectors are stored in
!     current reduced set).
!     -----------------------------------------------------------------

      IF (RUN_MODE .EQ. RUN_EXTERNAL) THEN
         CALL CHO_VECBUF_RETRIEVE(SCR,LSCR,JVEC1,IVEC2,ISYM,
     &                            JNUM,IREDC,MUSED)
      END IF

!     Read remaining vectors from disk.
!     ---------------------------------

      DOREAD = .TRUE.
      JV1 = JVEC1 + JNUM
      LFT = LSCR  - MUSED
      IF (IVEC2.GE.JV1 .AND. LFT.GT.0) THEN
         KS  = MUSED + 1
         JN  = 0
         MU  = 0
         CALL CHO_VECRD1(SCR(KS),LFT,JV1,IVEC2,ISYM,
     &                   JN,IREDC,MU,DOREAD)
         JNUM  = JNUM  + JN
         MUSED = MUSED + MU
      END IF

      END
      SUBROUTINE CHO_VECRD1(SCR,LSCR,JVEC1,IVEC2,ISYM,
     &                      JNUM,IREDC,MUSED,DOREAD)
!
!     Purpose: read as many vectors as fit into SCR array starting
!              at vector JVEC1 and reading at most until vector IVEC2.
!              On exit, JNUM is the number of vectors read.
!              On entry as well as exit, IREDC identifies the reduced
!              set stored in core (at position "3"; use -1 if none
!              or unkown). If DOREAD=.false. no vectors are actually
!              read in, but JNUM and MUSED are returned as appropriate.
!              Thus, array SCR is not referenced for DOREAD=.false.
!
!     NOTE: if no vectors can be read, JNUM=0 and MUSED=0 are returned,
!           but execution is NOT stopped here!!!
!
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
      Implicit Real*8 (a-h,o-z)
      DIMENSION SCR(LSCR)
      LOGICAL   DOREAD
#include "cholesky.fh"

      external ddot_

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_VECRD1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL FULL

      JRED = 0  ! fix compiler warning

      IF (CHO_ADRVEC .EQ. 1) THEN  ! WA addressing

!        Count how many vectors can be read.
!        -----------------------------------

         JNUM = 0
         LTOT = 0
         JVEC = JVEC1 - 1
         FULL = LTOT .GE. LSCR
         IF (.NOT.Allocated(NDIMRS)) THEN
            ILOC = 3 ! use scratch location in reduced index arrays
            DO WHILE ((JVEC.LT.IVEC2) .AND. (.NOT.FULL))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               IF (JRED .NE. IREDC) THEN
                  CALL CHO_GETRED(JRED,ILOC,.FALSE.)
                  CALL CHO_SETREDIND(ILOC)
                  IREDC = JRED
               END IF
               LTOT = LTOT + NNBSTR(ISYM,ILOC)
               IF (LTOT .GT. LSCR) THEN
                  JVEC = JVEC - 1
                  LTOT = LTOT - NNBSTR(ISYM,ILOC)
                  FULL = .TRUE.
               ELSE
                  JNUM = JNUM + 1
               END IF
            END DO
         ELSE
            DO WHILE ((JVEC.LT.IVEC2) .AND. (.NOT.FULL))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               LTOT = LTOT + NDIMRS(ISYM,JRED)
               IF (LTOT .GT. LSCR) THEN
                  JVEC = JVEC - 1
                  LTOT = LTOT - NDIMRS(ISYM,JRED)
                  FULL = .TRUE.
               ELSE
                  JNUM = JNUM + 1
               END IF
            END DO
         END IF

!        Read vectors (if any).
!        ----------------------

         IF (DOREAD .AND. LTOT.GT.0) THEN
            IOPT = 2
            IADR = INFVEC(JVEC1,3,ISYM)
            CALL DDAFILE(LUCHO(ISYM),IOPT,SCR,LTOT,IADR)
         END IF

      ELSE IF (CHO_ADRVEC .EQ. 2) THEN ! DA adressing

!        Read as many vectors as can be read, one at a time.
!        ---------------------------------------------------

         JNUM = 0
         LTOT = 0
         KSCR = 1
         JVEC = JVEC1 - 1
         FULL = LTOT .GE. LSCR
         IF (.NOT.Allocated(NDIMRS)) THEN
            ILOC = 3 ! use scratch location in reduced index arrays
            DO WHILE ((JVEC.LT.IVEC2) .AND. (.NOT.FULL))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               IF (JRED .NE. IREDC) THEN
                  CALL CHO_GETRED(JRED,ILOC,.FALSE.)
                  CALL CHO_SETREDIND(ILOC)
                  IREDC = JRED
               END IF
               LTOT = LTOT + NNBSTR(ISYM,ILOC)
               IF (LTOT .GT. LSCR) THEN
                  JVEC = JVEC - 1
                  LTOT = LTOT - NNBSTR(ISYM,ILOC)
                  FULL = .TRUE.
               ELSE
                  JNUM = JNUM + 1
                  IF (DOREAD) THEN
                     IOPT = 2
                     LENR = NNBSTR(ISYM,ILOC)
                     IADR = INFVEC(JVEC,3,ISYM)
                     CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LENR,IADR)
                     KSCR = KSCR + NNBSTR(ISYM,ILOC)
                  END IF
               END IF
            END DO
         ELSE
            DO WHILE ((JVEC.LT.IVEC2) .AND. (.NOT.FULL))
               JVEC = JVEC + 1
               JRED = INFVEC(JVEC,2,ISYM)
               LTOT = LTOT + NDIMRS(ISYM,JRED)
               IF (LTOT .GT. LSCR) THEN
                  JVEC = JVEC - 1
                  LTOT = LTOT - NDIMRS(ISYM,JRED)
                  FULL = .TRUE.
               ELSE
                  JNUM = JNUM + 1
                  IF (DOREAD) THEN
                     IOPT = 2
                     LENR = NDIMRS(ISYM,JRED)
                     IADR = INFVEC(JVEC,3,ISYM)
                     CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LENR,IADR)
                     KSCR = KSCR + NDIMRS(ISYM,JRED)
                  END IF
               END IF
            END DO
         END IF

      ELSE  ! unknown file addressing

         CALL CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
         LTOT = 0 ! dummy assignment to avoid compiler warnings

      END IF

!     Return total memory used for read.
!     ----------------------------------

      MUSED = LTOT

!     Debug print.
!     ------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Vector addressing: ',CHO_ADRVEC
         WRITE(LUPRI,*) 'DOREAD: ',DOREAD
         IF (JNUM .LT. 1) THEN
            IF (DOREAD) THEN
               WRITE(LUPRI,*) 'No vectors read!'
            ELSE
               WRITE(LUPRI,*) 'No vectors can be read!'
            END IF
         ELSE
            IF (DOREAD) THEN
               WRITE(LUPRI,*) 'Vectors ',JVEC1,' to ',JVEC1+JNUM-1,
     &                        ' of symmetry ',ISYM,' read from unit ',
     &                        LUCHO(ISYM)
               IF (Allocated(NDIMRS)) THEN
                  KOFFV = 1
                  DO IVEC = 1,JNUM
                     JVEC = JVEC1 + IVEC - 1
                     JADR = INFVEC(JVEC,3,ISYM)
                     JRED = INFVEC(JVEC,2,ISYM)
                     XNRM = SQRT(DDOT_(NDIMRS(ISYM,JRED),SCR(KOFFV),1,
     &                                                  SCR(KOFFV),1))
                     WRITE(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,
     &                              ' norm: ',XNRM
                     KOFFV = KOFFV + NDIMRS(ISYM,JRED)
                  END DO
                  NTST = KOFFV - 1
                  IF (NTST .NE. MUSED) THEN
                    CALL CHO_QUIT('Vector dimension error in '//SECNAM,
     &                            104)
                  END IF
               END IF
            ELSE
               WRITE(LUPRI,*) 'Vectors ',JVEC1,' to ',JVEC1+JNUM-1,
     &                        ' of symmetry ',ISYM,' can be read'
            END IF
         END IF
         CALL CHO_FLUSH(LUPRI)
      END IF

      END
