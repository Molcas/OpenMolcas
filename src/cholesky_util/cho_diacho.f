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
      SUBROUTINE CHO_DIACHO(DIAG,ISYM,WRK,LWRK)
!
!     Purpose: update (i.e. subtract contributions from vectors on disk)
!              of symmetry block ISYM of diagonal in red. set 1.
!              This emulates the actual procedure during decomposition.
!
      use ChoSwp, only: InfVec, IndRed
      Implicit Real*8 (a-h,o-z)
      Real*8 Diag(*), WRK(LWRK)
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_DIACHO')

      LOGICAL SCDIAG_SAVE

      PARAMETER (N2 = INFVEC_N2)
      PARAMETER (ZERO = 0.0D0)

!     Return if nothing to do.
!     ------------------------

      IF (NNBSTR(ISYM,1) .LT. 1) RETURN
      IF (NUMCHO(ISYM)   .LT. 1) RETURN

!     Save read-call counter.
!     -----------------------

      NSCALL = NSYS_CALL

!     Set pointer to scratch for reduced set indices.
!     -----------------------------------------------

      ILOC  = 3

!     Set up rs1 indices at location ILOC.
!     Set IREDC to identify this.
!     ------------------------------------

      CALL CHO_RSCOPY(1,ILOC)
      IREDC = 1

!     Start read buffer batch loop.
!     -----------------------------

      IVEC1 = 1
      DO WHILE (IVEC1 .LE. NUMCHO(ISYM))

!        Read as many vectors as possible into buffer (entire WRK).
!        ----------------------------------------------------------

         NVRD  = 0
         MUSED = 0
         CALL CHO_VECRD(WRK,LWRK,IVEC1,NUMCHO(ISYM),ISYM,
     &                  NVRD,IREDC,MUSED)
         IF (NVRD .LT. 1) THEN
            CALL CHO_QUIT('Insufficient scratch space for read in '
     &                    //SECNAM,101)
         END IF

!        Initialize vector offset.
!        -------------------------

         KOFFV = 0

!        Loop over vectors in core.
!        --------------------------

         DO JVEC = 1,NVRD

!           Set index arrays for current reduced set (if not already
!           set).
!           --------------------------------------------------------

            JRED = INFVEC(IVEC1+JVEC-1,2,ISYM)
            IF (JRED .NE. IREDC) THEN
               IF (JRED .EQ. 1) THEN
                  CALL CHO_RSCOPY(1,ILOC)
               ELSE
                  CALL CHO_GETRED(JRED,ILOC,.FALSE.)
                  CALL CHO_SETREDIND(ILOC)
               END IF
               IREDC = JRED
            END IF

!           Compute contributions to diagonal.
!           Zero the diagonal element associated with this vector.
!           ------------------------------------------------------

            DO JAB = 1,NNBSTR(ISYM,ILOC)
               IAB = INDRED(IIBSTR(ISYM,ILOC)+JAB,ILOC) ! address in rs1
               KAB = KOFFV + JAB ! vector address
               DIAG(IAB) = DIAG(IAB) - WRK(KAB)*WRK(KAB)
            END DO
            IABG = INFVEC(IVEC1+JVEC-1,1,ISYM)
            CALL CHO_P_ZERODIAG_RST(DIAG,ISYM,IABG)

!           Check diagonal.
!           ---------------

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


!           Update vector offset.
!           ---------------------

            KOFFV = KOFFV + NNBSTR(ISYM,ILOC)

         END DO

!        Check memory.
!        -------------

         IF (KOFFV .NE. MUSED) THEN
            CALL CHO_QUIT('Memory error detected in '//SECNAM,101)
         END IF

!        Update vector counter.
!        ----------------------

         IVEC1 = IVEC1 + NVRD

      END DO

!     Restore read-call counter.
!     --------------------------

      NSYS_CALL = NSCALL

      END
