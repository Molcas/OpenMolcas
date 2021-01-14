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
      SUBROUTINE CHO_DBGINT_CHO(XINT,NCD,NAB,WRK,LWRK,
     &                          ERRMAX,ERRMIN,ERRRMS,NCMP,
     &                          ISHLCD,ISHLAB)
C
C     Purpose: calculate integrals in shell quadruple (CD|AB) from
C              Cholesky vectors on disk and compare to those in
C              XINT (for debugging).
C
C     NOTE: this is *only* for debugging.
C
      use ChoArr, only: iSP2F
#include "implicit.fh"
      DIMENSION XINT(NCD,NAB), WRK(LWRK)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DBGINT_CHO')

      INTEGER  CHO_LREAD
      EXTERNAL CHO_LREAD

      IIBSTRSH(I,J,K)=IWORK(ip_IIBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      INDRED(I,J)=IWORK(ip_INDRED-1+MMBSTRT*(J-1)+I)
      NBSTSH(I)=IWORK(ip_NBSTSH-1+I)

C     Initializations.
C     ----------------


      CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
      ERRMAX = -1.0D12
      ERRMIN =  1.0D12
      ERRRMS = 0.0D0
      NCMP   = 0
      LCDABT = NCD*NAB

      IF (ISHLC .EQ. ISHLD) THEN
         NCDL = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
      ELSE
         NCDL = NBSTSH(ISHLC)*NBSTSH(ISHLD)
      END IF
      IF (ISHLA .EQ. ISHLB) THEN
         NABL = NBSTSH(ISHLA)*(NBSTSH(ISHLA) + 1)/2
      ELSE
         NABL = NBSTSH(ISHLA)*NBSTSH(ISHLB)
      END IF
      IF (NCDL .GT. NCD) CALL CHO_QUIT('NCD error in '//SECNAM,104)
      IF (NABL .GT. NAB) CALL CHO_QUIT('NAB error in '//SECNAM,104)
      IF (NAB.LT.1 .OR. NCD.LT.1) RETURN

C     Save read-call counter.
C     -----------------------

      NSCALL = NSYS_CALL

C     Get a copy of XINT.
C     -------------------

      KXINT = 1
      KEND0 = KXINT + LCDABT
      LWRK0 = LWRK  - KEND0
      IF (LWRK0 .LE. 0) THEN
         CALL CHO_QUIT('Insufficient memory in '//SECNAM//' [0]',101)
      END IF

      CALL DCOPY_(LCDABT,XINT,1,WRK(KXINT),1)

C     Start symmetry loop.
C     --------------------

      DO ISYM = 1,NSYM

         NUMCD = NNBSTRSH(ISYM,ISHLCD,2)
         NUMAB = NNBSTRSH(ISYM,ISHLAB,2)

         IF (NUMCD.GT.0 .AND. NUMAB.GT.0 .AND. NUMCHO(ISYM).GT.0) THEN

C           Allocate space for integrals and for Cholesky reading.
C           ------------------------------------------------------

            LENint = NUMCD*NUMAB
            LREAD  = CHO_LREAD(ISYM,LWRK)
            LVEC1  = NNBSTR(ISYM,2)

            KINT  = KEND0
            KREAD = KINT  + LENint
            KVEC1 = KREAD + LREAD
            KEND1 = KVEC1 + LVEC1
            LWRK1 = LWRK  - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               CALL CHO_QUIT('Insufficient memory in '//SECNAM,104)
            END IF

C           Initialize integral array.
C           --------------------------

            CALL CHO_DZERO(WRK(KINT),LENint)

C           Set up batch over Cholesky vectors.
C           -----------------------------------

            MINM = NUMCD + NUMAB
            NVEC = MIN(LWRK1/MINM,NUMCHO(ISYM))
            IF (NVEC .LT. 1) THEN
               CALL CHO_QUIT('Batch problem in '//SECNAM,104)
            END IF
            NBATCH = (NUMCHO(ISYM) - 1)/NVEC + 1

C           Start batch loop.
C           -----------------

            DO IBATCH = 1,NBATCH

               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NUMCHO(ISYM) - NVEC*(NBATCH - 1)
               ELSE
                  NUMV = NVEC
               END IF
               JVEC1 = NVEC*(IBATCH - 1) + 1

               KCHOCD = KEND1
               KCHOAB = KCHOCD + NUMCD*NUMV
               KEND2  = KCHOAB + NUMAB*NUMV
               LWRK2  = LWRK   - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL CHO_QUIT('Batch error in '//SECNAM,104)
               END IF

C              Read vectors.
C              -------------

               DO IVEC = 1,NUMV
                  JVEC = JVEC1 + IVEC - 1
                  CALL CHO_GETVEC(WRK(KVEC1),LVEC1,1,JVEC,ISYM,
     &                            WRK(KREAD),LREAD)
                  KOFF1 = KVEC1  + IIBSTRSH(ISYM,ISHLCD,2)
                  KOFF2 = KCHOCD + NUMCD*(IVEC - 1)
                  CALL DCOPY_(NUMCD,WRK(KOFF1),1,WRK(KOFF2),1)
                  KOFF1 = KVEC1  + IIBSTRSH(ISYM,ISHLAB,2)
                  KOFF2 = KCHOAB + NUMAB*(IVEC - 1)
                  CALL DCOPY_(NUMAB,WRK(KOFF1),1,WRK(KOFF2),1)
               END DO

C              Calculate contribution.
C              -----------------------

               CALL DGEMM_('N','T',NUMCD,NUMAB,NUMV,
     &                    1.0D0,WRK(KCHOCD),NUMCD,WRK(KCHOAB),NUMAB,
     &                    1.0D0,WRK(KINT),NUMCD)

            END DO

C           Subtract contribution from full shell pair.
C           -------------------------------------------

            DO IAB = 1,NUMAB
               JAB = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLAB,2) + IAB
               KAB = INDRED(INDRED(JAB,2),1)
               DO ICD = 1,NUMCD
                  JCD = IIBSTR(ISYM,2) + IIBSTRSH(ISYM,ISHLCD,2) + ICD
                  KCD = INDRED(INDRED(JCD,2),1)
                  ICDAB = KINT  + NUMCD*(IAB - 1) + ICD - 1
                  KCDAB = KXINT + NCD*(KAB - 1)   + KCD - 1
                  WRK(KCDAB) = WRK(KCDAB) - WRK(ICDAB)
               END DO
            END DO

         END IF

      END DO

C     Compare full shell pair.
C     ------------------------

      DO KAB = 1,NAB
         DO KCD = 1,NCD
            KCDAB = KXINT + NCD*(KAB - 1) + KCD - 1
            DIFF  = WRK(KCDAB)
            NCMP  = NCMP + 1
            IF (NCMP .EQ. 1) THEN
               ERRMAX = DIFF
               ERRMIN = DIFF
            ELSE
               IF (ABS(DIFF) .GT. ABS(ERRMAX)) THEN
                  ERRMAX = DIFF
               END IF
               IF (ABS(DIFF) .LT. ABS(ERRMIN)) THEN
                  ERRMIN = DIFF
               END IF
            END IF
            ERRRMS = ERRRMS + DIFF*DIFF
         END DO
      END DO

C     Restore read-call counter.
C     --------------------------

      NSYS_CALL = NSCALL

      END
