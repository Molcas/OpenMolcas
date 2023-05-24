!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************
      SUBROUTINE CHO_SUBTR0(XINT,WRK,LWRK,ISYM)
!
!     Purpose: subtract contributions from previous vectors
!              from the qualified integrals (in XINT).
!              This version is memory-driven.
!
!     Screening in subtraction introduced Jan. 2006, TBP.
!
      use ChoSwp, only: iQuAB, nnBstRSh, iiBstRSh
      use ChoArr, only: LQ
      use ChoVecBuf, only: nVec_in_Buf
      use ChoSubScr, only: Cho_SScreen, SSTau, SubScrStat, DSubScr,
     &                     DSPNm, SSNorm
      Implicit Real*8 (a-h,o-z)
      DIMENSION XINT(*), WRK(LWRK)
#include "cholesky.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_SUBTR0')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER  CHO_LREAD
      EXTERNAL CHO_LREAD

      PARAMETER (XMONE = -1.0D0, ONE = 1.0D0)

!     Return if nothing to do.
!     ------------------------

      IF (NUMCHO(ISYM) .LT. 1) RETURN

      NVEC_TO_READ = NUMCHO(ISYM) - NVEC_IN_BUF(ISYM)
      IF (NVEC_TO_READ .EQ. 0) RETURN
      IF (NVEC_TO_READ .LT. 0) THEN
         CALL CHO_QUIT('Vector buffer error in '//SECNAM,104)
      END IF

!     Initialize.
!     -----------

      XTOT = 0.0D0
      XDON = 0.0D0

!     Reserve space needed for reading previous vectors.
!     --------------------------------------------------

      LREAD = CHO_LREAD(ISYM,LWRK)
      IF (LREAD .LT. 1) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_LREAD returned ',LREAD
         CALL CHO_QUIT('Memory error in '//SECNAM,101)
         LWRK1 = 0 ! to avoid compiler warnings
      ELSE
         LWRK1 = LWRK - LREAD
      END IF

!     Set up batch.
!     -------------

      MMEM = NNBSTR(ISYM,2) + NQUAL(ISYM)
      NVEC = MIN(LWRK1/MMEM,NVEC_TO_READ)
      IF (NVEC .LT. 1) THEN
         CALL CHO_QUIT('Batch failure in '//SECNAM,101)
      END IF
      NBATCH = (NVEC_TO_READ - 1)/NVEC + 1

!     Start batch loop.
!     -----------------

      DO IBATCH = 1,NBATCH

         IF (IBATCH .EQ. NBATCH) THEN
            NUMV = NVEC_TO_READ - NVEC*(NBATCH - 1)
         ELSE
            NUMV = NVEC
         END IF
         IVEC1 = NVEC_IN_BUF(ISYM) + NVEC*(IBATCH - 1) + 1

!        Complete allocation.
!        --------------------

         KCHO1 = 1
         KCHO2 = KCHO1 + NNBSTR(ISYM,2)*NUMV
         KEND2 = KCHO2 + NQUAL(ISYM)*NUMV
         LWRK2 = LWRK  - KEND2 + 1
         IF (LWRK2 .LT. LREAD) THEN
            CALL CHO_QUIT('Batch error in '//SECNAM,104)
         END IF

!        Read previous vectors.
!        ----------------------

         CALL CHO_TIMER(C1,W1)
         CALL CHO_GETVEC(WRK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,ISYM,
     &                   WRK(KEND2),LWRK2)
         CALL CHO_TIMER(C2,W2)
         TDECOM(1,2) = TDECOM(1,2) + C2 - C1
         TDECOM(2,2) = TDECOM(2,2) + W2 - W1


!        Screened or unscreened subtraction section.
!        The screened version uses level 2 blas, while the unscreened
!        one employs level 3 blas.
!        ------------------------------------------------------------

         CALL CHO_TIMER(C1,W1)

         IF (CHO_SSCREEN) THEN ! screened subtraction

            KOFB0 = KCHO1 - 1 - IIBSTR(ISYM,2)
            DO J = 1,NUMV
               KOFFA = KCHO2 + J - 1
               KOFFB = KOFB0 + NNBSTR(ISYM,2)*(J - 1)
               DO IAB = 1,NQUAL(ISYM)
                  WRK(KOFFA+NUMV*(IAB-1))=WRK(KOFFB+IQUAB(IAB,ISYM))
               END DO
            END DO

!           Subtract:
!           (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L(#J,{ab})
!           for each ab in {ab}.
!           ----------------------------------------------------

            CALL CHO_SUBSCR_DIA(WRK(KCHO1),NUMV,iSym,2,SSNorm)
            DO IAB = 1,NQUAL(ISYM)
               DO ISHGD = 1,NNSHL
                  NGD = NNBSTRSH(ISYM,ISHGD,2)
                  IF (NGD .GT. 0) THEN
                     XTOT = XTOT + 1.0D0
                     JAB = IQUAB(IAB,ISYM) - IIBSTR(ISYM,2)
                     TST = SQRT(DSPNM(ISHGD)*DSUBSCR(JAB))
                     IF (TST .GT. SSTAU) THEN
                        XDON = XDON + 1.0D0
                        KOFF1 = KCHO1 + IIBSTRSH(ISYM,ISHGD,2)
                        KOFF2 = KCHO2 + NUMV*(IAB-1)
                        KOFF3 = NNBSTR(ISYM,2)*(IAB-1)
     &                        + IIBSTRSH(ISYM,ISHGD,2) + 1
                        CALL DGEMV_('N',NGD,NUMV,
     &                             XMONE,WRK(KOFF1),NNBSTR(ISYM,2),
     &                             WRK(KOFF2),1,ONE,XINT(KOFF3),1)
                     END IF
                  END IF
               END DO
            END DO

         ELSE ! unscreened subtraction

           IF (Associated(LQ(ISYM)%Array)) THEN

!              If the qualified block, L({ab},#J), is already in core,
!              use this block.
!              -------------------------------------------------------

               CALL DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,
     &                    XMONE,WRK(KCHO1),NNBSTR(ISYM,2),
     &                          LQ(ISYM)%Array(:,IVEC1),
     &                          SIZE(LQ(ISYM)%Array,1),
     &                    ONE,XINT,NNBSTR(ISYM,2))


            ELSE

!              Copy out sub-blocks corresponding to qualified diagonals:
!              L({ab},#J)
!              ---------------------------------------------------------

               DO J = 1,NUMV
                  DO IAB = 1,NQUAL(ISYM)
                     KOFF1 = KCHO2 + NQUAL(ISYM)*(J - 1)
     &                     + IAB - 1
                     KOFF2 = KCHO1 + NNBSTR(ISYM,2)*(J - 1)
     &                     + IQUAB(IAB,ISYM) - IIBSTR(ISYM,2) - 1
                     WRK(KOFF1) = WRK(KOFF2)
                  END DO
               END DO

!              Subtract:
!              (gd|{ab}) <- (gd|{ab}) - sum_J L(gd,#J) * L({ab},#J)
!              ----------------------------------------------------

               CALL DGEMM_('N','T',NNBSTR(ISYM,2),NQUAL(ISYM),NUMV,
     &                    XMONE,WRK(KCHO1),NNBSTR(ISYM,2),
     &                          WRK(KCHO2),NQUAL(ISYM),
     &                    ONE,XINT,NNBSTR(ISYM,2))

            END IF

!           Update DGEMM-call counter.
!           --------------------------

            NDGM_CALL = NDGM_CALL + 1

         END IF

         CALL CHO_TIMER(C2,W2)
         TDECOM(1,3) = TDECOM(1,3) + C2 - C1
         TDECOM(2,3) = TDECOM(2,3) + W2 - W1

      END DO

!     Update screening statistics.
!     ----------------------------

      IF (CHO_SSCREEN) THEN
         SUBSCRSTAT(1) = SUBSCRSTAT(1) + XTOT
         SUBSCRSTAT(2) = SUBSCRSTAT(2) + XDON
      END IF

      END
