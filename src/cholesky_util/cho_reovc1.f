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
      SUBROUTINE CHO_REOVC1(IRS2F,N,LRDIM,WRK,LWRK)
C
C     Purpose: reorder Cholesky vectors on disk to full storage.
C
      use ChoSwp, only: nnBstRSh
#include "implicit.fh"
      INTEGER   IRS2F(N,LRDIM)
      DIMENSION WRK(LWRK)
#include "cholesky.fh"
#include "choreo.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_REOVC1')

      INTEGER IOFF(8,8)

      MULD2H(I,J)=IEOR(I-1,J-1)+1

      IF (N .LT. 3) THEN
         CALL CHO_QUIT('Dimension error in '//SECNAM,104)
      END IF

C     Save read-call counter.
C     -----------------------

      NSCALL = NSYS_CALL

C     Make rs1 the "current" reduced set (for reading).
C     -------------------------------------------------

      CALL CHO_RSCOPY(IWORK(ip_IIBSTRSH),NNBSTRSH,
     &                IWORK(ip_INDRED),1,2,NSYM,NNSHL,NNBSTRT(1),3)

C     Loop over Cholesky vector symmetries.
C     -------------------------------------

      DO ISYM = 1,NSYM
         IF (NUMCHO(ISYM) .GT. 0) THEN

C           Open files.
C           -----------

            CALL CHO_OPFVEC(ISYM,1)

C           Set up vector batch.
C           --------------------

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

C           Start batch loop over vectors.
C           ------------------------------

            DO IBATCH = 1,NBATCH

               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
               ELSE
                  NUMV = NVEC
               END IF
               IVEC1 = NVEC*(IBATCH - 1) + 1

C              Read batch of reduced vectors.
C              ------------------------------

               KCHO1 = 1
               KREAD = KCHO1 + NNBSTR(ISYM,2)*NUMV
               LREAD = LWRK - KREAD + 1
               CALL CHO_GETVEC(WRK(KCHO1),NNBSTR(ISYM,2),NUMV,IVEC1,
     &                         ISYM,WRK(KREAD),LREAD)

C              Reorder.
C              --------

               KCHO2 = KREAD
               CALL CHO_IZERO(IOFF,64)
               ICOUNT = KCHO2 - 1
               DO ISYMB = 1,NSYM
                  ISYMA = MULD2H(ISYMB,ISYM)
                  IF (ISYMA .GE. ISYMB) THEN
                     IOFF(ISYMA,ISYMB) = ICOUNT
                     IOFF(ISYMB,ISYMA) = ICOUNT
                     ICOUNT = ICOUNT + NABPK(ISYMA,ISYMB)*NUMV
                  END IF
               END DO

               CALL CHO_DZERO(WRK(KCHO2),NNBST(ISYM)*NUMV)
               DO IVEC = 1,NUMV
                  KOFF1 = KCHO1 + NNBSTR(ISYM,2)*(IVEC - 1) - 1
                  DO IRS = 1,NNBSTR(ISYM,2)
                     I = IIBSTR(ISYM,2) + IRS
                     ISYMA = IRS2F(1,I)
                     ISYMB = IRS2F(2,I)
                     IAB   = IRS2F(3,I)
                     KOFF  = KOFF1 + IRS
                     LOFF  = IOFF(ISYMA,ISYMB)
     &                     + NABPK(ISYMA,ISYMB)*(IVEC - 1) + IAB
                     WRK(LOFF) = WRK(KOFF)
                  END DO
               END DO

C              Write full vectors to disk.
C              ---------------------------

               DO ISYMB = 1,NSYM
                  ISYMA = MULD2H(ISYMB,ISYM)
                  IF (ISYMA .GE. ISYMB) THEN
                     KOFF = IOFF(ISYMA,ISYMB) + 1
                     CALL CHO_WRFVEC(WRK(KOFF),ISYMA,ISYMB,IVEC1,NUMV)
                  END IF
               END DO

            END DO

C           Close files.
C           ------------

            CALL CHO_OPFVEC(ISYM,2)

         END IF
      END DO

C     Restore read-call counter.
C     --------------------------

      NSYS_CALL = NSCALL

      END
