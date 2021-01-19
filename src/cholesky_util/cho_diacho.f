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
      SUBROUTINE CHO_DIACHO(DIAG,ISYM,WRK,LWRK)
C
C     Purpose: update (i.e. subtract contributions from vectors on disk)
C              of symmetry block ISYM of diagonal in red. set 1.
C              This emulates the actual procedure during decomposition.
C
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, InfRed, InfVec,
     &                  IndRed
#include "implicit.fh"
      DIMENSION DIAG(*), WRK(LWRK)
#include "cholesky.fh"
#include "choptr.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_DIACHO')

      LOGICAL SCDIAG_SAVE

      PARAMETER (N2 = INFVEC_N2)
      PARAMETER (ZERO = 0.0D0)

C     Return if nothing to do.
C     ------------------------

      IF (NNBSTR(ISYM,1) .LT. 1) RETURN
      IF (NUMCHO(ISYM)   .LT. 1) RETURN

C     Save read-call counter.
C     -----------------------

      NSCALL = NSYS_CALL

C     Set pointer to scratch for reduced set indices.
C     -----------------------------------------------

      ILOC  = 3

C     Set up rs1 indices at location ILOC.
C     Set IREDC to identify this.
C     ------------------------------------

      CALL CHO_RSCOPY(IIBSTRSH,NNBSTRSH,
     &                INDRED,1,ILOC,NSYM,NNSHL,NNBSTRT(1),3)
      IREDC = 1

C     Start read buffer batch loop.
C     -----------------------------

      IVEC1 = 1
      DO WHILE (IVEC1 .LE. NUMCHO(ISYM))

C        Read as many vectors as possible into buffer (entire WRK).
C        ----------------------------------------------------------

         NVRD  = 0
         MUSED = 0
         CALL CHO_VECRD(WRK,LWRK,IVEC1,NUMCHO(ISYM),ISYM,
     &                  NVRD,IREDC,MUSED)
         IF (NVRD .LT. 1) THEN
            CALL CHO_QUIT('Insufficient scratch space for read in '
     &                    //SECNAM,101)
         END IF

C        Initialize vector offset.
C        -------------------------

         KOFFV = 0

C        Loop over vectors in core.
C        --------------------------

         DO JVEC = 1,NVRD

C           Set index arrays for current reduced set (if not already
C           set).
C           --------------------------------------------------------

            JRED = INFVEC(IVEC1+JVEC-1,2,ISYM)
            IF (JRED .NE. IREDC) THEN
               IF (JRED .EQ. 1) THEN
                  CALL CHO_RSCOPY(IIBSTRSH,NNBSTRSH,INDRED,1,ILOC,
     &                            NSYM,NNSHL,NNBSTRT(1),3)
               ELSE
                  CALL CHO_GETRED(INFRED,nnBstRSh(:,:,ILOC),
     &                            IndRed(1,ILOC),INDRSH,iSP2F,
     &                            MAXRED,NSYM,NNSHL,MMBSTRT,JRED,
     &                            .FALSE.)
                  CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,ILOC)
               END IF
               IREDC = JRED
            END IF

C           Compute contributions to diagonal.
C           Zero the diagonal element associated with this vector.
C           ------------------------------------------------------

            DO JAB = 1,NNBSTR(ISYM,ILOC)
               IAB = INDRED(IIBSTR(ISYM,ILOC)+JAB,ILOC) ! address in rs1
               KAB = KOFFV + JAB ! vector address
               DIAG(IAB) = DIAG(IAB) - WRK(KAB)*WRK(KAB)
            END DO
            IABG = INFVEC(IVEC1+JVEC-1,1,ISYM)
            CALL CHO_P_ZERODIAG_RST(DIAG,ISYM,IABG)

C           Check diagonal.
C           ---------------

            IF (CHO_DECALG .EQ. 4) THEN
               SCDIAG_SAVE = SCDIAG
               SCDIAG = .FALSE. ! do NOT screen
               DMX = 1.0D0
               CALL CHO_CHKDIA_A4(DIAG,DMX,ISYM,NNEG,NNEGT,NCONV,XMAX,
     &                            XMIN,XM)
               SCDIAG = SCDIAG_SAVE
            ELSE
               CALL CHO_CHKDIA(DIAG,ISYM,XMIN,XMAX,XM,NNEGT,NNEG,NCONV)
            END IF


C           Update vector offset.
C           ---------------------

            KOFFV = KOFFV + NNBSTR(ISYM,ILOC)

         END DO

C        Check memory.
C        -------------

         IF (KOFFV .NE. MUSED) THEN
            CALL CHO_QUIT('Memory error detected in '//SECNAM,101)
         END IF

C        Update vector counter.
C        ----------------------

         IVEC1 = IVEC1 + NVRD

      END DO

C     Restore read-call counter.
C     --------------------------

      NSYS_CALL = NSCALL

      END
