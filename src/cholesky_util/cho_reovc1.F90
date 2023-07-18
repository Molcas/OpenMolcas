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
      SUBROUTINE CHO_REOVC1(IRS2F,N,LRDIM,WRK,LWRK)
      use ChoReO
!
!     Purpose: reorder Cholesky vectors on disk to full storage.
!
      Implicit Real*8 (a-h,o-z)
      INTEGER   N,LRDIM,LWRK
      INTEGER   IRS2F(N,LRDIM)
      REAL*8 WRK(LWRK)
#include "cholesky.fh"

      CHARACTER(LEN=10), PARAMETER::SECNAM = 'CHO_REOVC1'
      INTEGER IOFF(8,8), I, J, MulD2h

      MULD2H(I,J)=IEOR(I-1,J-1)+1

      IF (N .LT. 3) THEN
         CALL CHO_QUIT('Dimension error in '//SECNAM,104)
      END IF

!     Save read-call counter.
!     -----------------------

      NSCALL = NSYS_CALL

!     Make rs1 the "current" reduced set (for reading).
!     -------------------------------------------------

      CALL CHO_RSCOPY(1,2)

!     Loop over Cholesky vector symmetries.
!     -------------------------------------

      DO ISYM = 1,NSYM
         IF (NUMCHO(ISYM) .GT. 0) THEN

!           Open files.
!           -----------

            CALL CHO_OPFVEC(ISYM,1)

!           Set up vector batch.
!           --------------------

            MINMEM = NNBSTR(ISYM,2) + NNBST(ISYM)
            IF (MINMEM .LT. 1) THEN
               WRITE(LUPRI,*) SECNAM,': MINMEM = ',MINMEM
               CALL CHO_QUIT('NNBST error in '//SECNAM,104)
               NVEC = 0
            ELSE
               NVEC = MIN(LWRK/MINMEM,NUMCHO(ISYM))
            END IF

            IF (NVEC .LT. 1) THEN
               WRITE(LUPRI,*) SECNAM,': NVEC   = ',NVEC
               WRITE(LUPRI,*) SECNAM,': LWRK   = ',LWRK
               WRITE(LUPRI,*) SECNAM,': MINMEM = ',MINMEM
               WRITE(LUPRI,*) SECNAM,': NUMCHO = ',NUMCHO(ISYM)
               WRITE(LUPRI,*) SECNAM,': ISYM   = ',ISYM
               CALL CHO_QUIT('Batch error in '//SECNAM,101)
               NBATCH = 0
            ELSE
               NBATCH = (NUMCHO(ISYM) - 1)/NVEC + 1
            END IF

!           Start batch loop over vectors.
!           ------------------------------

            DO IBATCH = 1,NBATCH

               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
               ELSE
                  NUMV = NVEC
               END IF
               IVEC1 = NVEC*(IBATCH - 1) + 1

!              Read batch of reduced vectors.
!              ------------------------------

               KCHO1 = 1
               KREAD = KCHO1 + NNBSTR(ISYM,2)*NUMV
               LREAD = LWRK - KREAD + 1
               CALL CHO_GETVEC(WRK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,    &
     &                         ISYM,WRK(KREAD),LREAD)

!              Reorder.
!              --------

               KCHO2 = KREAD
               CALL IZERO(IOFF,64)
               ICOUNT = KCHO2 - 1
               DO ISYMB = 1,NSYM
                  ISYMA = MULD2H(ISYMB,ISYM)
                  IF (ISYMA .GE. ISYMB) THEN
                     IOFF(ISYMA,ISYMB) = ICOUNT
                     IOFF(ISYMB,ISYMA) = ICOUNT
                     ICOUNT = ICOUNT + NABPK(ISYMA,ISYMB)*NUMV
                  END IF
               END DO

               CALL FZERO(WRK(KCHO2),NNBST(ISYM)*NUMV)
               DO IVEC = 1,NUMV
                  KOFF1 = KCHO1 + NNBSTR(ISYM,2)*(IVEC - 1) - 1
                  DO IRS = 1,NNBSTR(ISYM,2)
                     I = IIBSTR(ISYM,2) + IRS
                     ISYMA = IRS2F(1,I)
                     ISYMB = IRS2F(2,I)
                     IAB   = IRS2F(3,I)
                     KOFF  = KOFF1 + IRS
                     LOFF  = IOFF(ISYMA,ISYMB)                          &
     &                     + NABPK(ISYMA,ISYMB)*(IVEC - 1) + IAB
                     WRK(LOFF) = WRK(KOFF)
                  END DO
               END DO

!              Write full vectors to disk.
!              ---------------------------

               DO ISYMB = 1,NSYM
                  ISYMA = MULD2H(ISYMB,ISYM)
                  IF (ISYMA .GE. ISYMB) THEN
                     KOFF = IOFF(ISYMA,ISYMB) + 1
                     CALL CHO_WRFVEC(WRK(KOFF),ISYMA,ISYMB,IVEC1,NUMV)
                  END IF
               END DO

            END DO

!           Close files.
!           ------------

            CALL CHO_OPFVEC(ISYM,2)

         END IF
      END DO

!     Restore read-call counter.
!     --------------------------

      NSYS_CALL = NSCALL

      END
