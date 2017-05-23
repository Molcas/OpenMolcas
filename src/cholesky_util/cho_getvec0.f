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
      SUBROUTINE CHO_GETVEC0(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,
     &                       SCR,LSCR)
C
C=======================================================================
C==== DEPRECATED - USE CHO_X_GETVEC OR CHO_X_VECRD INSTEAD =============
C=======================================================================
C
C     Purpose: read Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
C              of symmetry ISYM from file. The vectors are returned
C              in the "current" reduced set.
C
C     NOTE: array SCR(LSCR) is used for storing the vectors in the
C           red. set from disk and for a full first red. set vector.
C           Thus, to be certain that enough memory is available,
C           use LSCR = 2 x dimension of first reduced set.
C
C
#include "implicit.fh"
      DIMENSION CHOVEC(LENVEC,NUMVEC)
      DIMENSION SCR(LSCR)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      external ddot_

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVEC0')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (N2 = INFVEC_N2)

      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)

C     Initialize output array.
C     ------------------------

      CALL CHO_DZERO(CHOVEC,LENVEC*NUMVEC)

C     Read reduced set index arrays for first vector.
C     -----------------------------------------------

      IRED  = INFVEC(IVEC1,2,ISYM)
      ILOC  = 3
      KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(ILOC - 1)
      KOFF2 = ip_INDRED   + MMBSTRT*(ILOC - 1)
      CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                IWORK(KOFF2),IWORK(ip_INDRSH),IWORK(ip_iSP2F),
     &                MAXRED,NSYM,NNSHL,MMBSTRT,IRED,
     &                .FALSE.)
      CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),
     &                   IWORK(ip_NNBSTRSH),NSYM,NNSHL,3)
      KRED1 = 1
      KREAD = KRED1 + NNBSTR(ISYM,1)
      KEND1 = KREAD + NNBSTR(ISYM,3)
      LSCR1 = LSCR  - KEND1 + 1
      IF (LSCR1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Insufficient scratch space in ',SECNAM
         WRITE(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
         WRITE(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
         WRITE(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
         CALL CHO_QUIT('[1] Insufficient scratch space in '//SECNAM,102)
      END IF

C     Read vectors and re-order into current reduced set via reduced
C     set 1.
C     NOTE: if the read vectors are already in red. set 1, don't resort.
C     ------------------------------------------------------------------

      DO JVEC = 1,NUMVEC

         IVEC = IVEC1 + JVEC - 1
         JRED = INFVEC(IVEC,2,ISYM)
         IF (JRED .NE. IRED) THEN   ! read new reduced set
            KOFF1 = ip_NNBSTRSH + NSYM*NNSHL*(ILOC - 1)
            KOFF2 = ip_INDRED   + MMBSTRT*(ILOC - 1)
            CALL CHO_GETRED(IWORK(ip_INFRED),IWORK(KOFF1),
     &                      IWORK(KOFF2),IWORK(ip_INDRSH),
     &                      IWORK(ip_iSP2F),
     &                      MAXRED,NSYM,NNSHL,MMBSTRT,JRED,
     &                      .FALSE.)
            CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),
     &                         IWORK(ip_NNBSTRSH),NSYM,NNSHL,3)
            KEND1 = KREAD + NNBSTR(ISYM,3)
            LSCR1 = LSCR  - KEND1 + 1
            IF (LSCR1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Insufficient scratch space in ',SECNAM
               WRITE(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
               WRITE(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
               WRITE(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
               CALL CHO_QUIT('[2] Insufficient scratch space in '
     &                       //SECNAM,102)
            END IF
            IRED = JRED
         END IF

C        IOPT = 2
C        IADR = INFVEC(IVEC,3,ISYM)
C        CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KREAD),NNBSTR(ISYM,3),IADR)
C-tbp: replaced above code to make use of buffer through cho_vecrd.
         JNUM  = 0
         IREDC = IRED
         MUSED = 0
         CALL CHO_VECRD(SCR(KREAD),NNBSTR(ISYM,3),IVEC,IVEC,ISYM,
     &                  JNUM,IREDC,MUSED)
         IF (JNUM .NE. 1) THEN
            CALL CHO_QUIT('Logical error in '//SECNAM,103)
         END IF
         NSYS_CALL = NSYS_CALL + 1
         IF (LOCDBG) THEN
            IADR = INFVEC(IVEC,3,ISYM)
            XNRM = SQRT(DDOT_(NNBSTR(ISYM,3),SCR(KREAD),1,SCR(KREAD),1))
            WRITE(LUPRI,*) SECNAM,': ',
     &                     'Vector:',IVEC,' address: ',
     &                     INFVEC(IVEC,3,ISYM),
     &                     ' norm: ',XNRM,' sym. ',ISYM,' red. set: ',
     &                     IRED,' dim.: ',NNBSTR(ISYM,3)
         END IF

         CALL CHO_DZERO(SCR(KRED1),NNBSTR(ISYM,1))
         IF (IRED .GT. 1) THEN
            DO JAB = 1,NNBSTR(ISYM,3)   ! sort into rs1 ordering
               KAB = IIBSTR(ISYM,3)  + JAB
               IAB = INDRED(KAB,3) - IIBSTR(ISYM,1)
               SCR(KRED1+IAB-1) = SCR(KREAD+JAB-1)
            END DO
            KREDU = KRED1   ! point rs2 sort to red1 resort
         ELSE IF (IRED .EQ. 1) THEN
            KREDU = KREAD   ! point rs2 sort to read (already red1)
         ELSE
            WRITE(LUPRI,*) SECNAM,': ERROR: IRED is negative: ',IRED
            CALL CHO_QUIT('Reduced set error in '//SECNAM,104)
            KREDU = -999999 ! just to avoid compiler warnings
         END IF

         DO JAB = 1,NNBSTR(ISYM,2)  ! sort into in rs2 ordering
            KAB = IIBSTR(ISYM,2)  + JAB
            IAB = INDRED(KAB,2) - IIBSTR(ISYM,1)
            CHOVEC(JAB,JVEC) = SCR(KREDU+IAB-1)
         END DO

      END DO

      END
