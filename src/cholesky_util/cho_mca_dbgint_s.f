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
      SUBROUTINE CHO_MCA_DBGINT_S(ISHLQ,NSHLQ,PRTLAB)
C
C     Purpose: regenerate and check a specified set of integrals.
C
C     NOTE: this is *not* meant for production calculations, only for
C           debugging, as:
C           1) full Cholesky vectors are read
C           2) calculations are performed in full (no use of red. sets
C              apart from first)
C
      use ChoArr, only: nBstSh
#include "implicit.fh"
      INTEGER ISHLQ(4,NSHLQ)
      LOGICAL PRTLAB
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      CHARACTER*16 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_DBGINT_S')

      INTEGER  CHO_F2SP
      EXTERNAL CHO_F2SP

      DIMENSION XLBAS(8)

      CHARACTER*8 LABEL

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      ITRI(I,J)=MAX(I,J)*(MAX(I,J)-3)/2+I+J

C     Return if nothing specified.
C     ----------------------------

      IF (NSHLQ .LT. 1) RETURN

C     Force computation of full shell quadruple.
C     ------------------------------------------

      IF (IFCSEW .NE. 1) THEN
         WRITE(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',
     &                  IFCSEW,' to 1.'
         WRITE(LUPRI,*) SECNAM,
     &   ': memory demands are significantly increased by this!'
         IFCSEW = 1
      END IF

C     Initializations.
C     ----------------

      GLMAX = 0.0D0
      GLMIN = 1.0D15
      GLRMS = 0.0D0
      XTCMP = 0.0D0
      XPECT = 0.0D0

C     Make first reduced set the current reduced set.
C     -----------------------------------------------

      CALL CHO_RSCOPY(1,2)

C     Allocate memory for largest integral quadruple.
C     -----------------------------------------------

      LINTMX = MX2SH*MX2SH
      CALL GETMEM('DBGINT.1','ALLO','REAL',KINT1,LINTMX)

C     Allocate max. memory
C     ----------------------------------------------------------

      CALL GETMEM('DBGINT.2','MAX ','REAL',KWRK,LWRK)
      CALL GETMEM('DBGINT.2','ALLO','REAL',KWRK,LWRK/2)
      CALL XSETMEM_INTS(LWRK/2)

C     Print header.
C     -------------

      CALL CHO_HEAD('Integral Error Analysis','=',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A)')
     & '    C     D     A     B   Abs. Min.    Abs. Max.      RMS',
     & '--------------------------------------------------------------'

C     Loop over specified shell quadruples.
C     -------------------------------------

      DO I = 1,NSHLQ

         ISHLC = ISHLQ(1,I)
         ISHLD = ISHLQ(2,I)
         ISHLA = ISHLQ(3,I)
         ISHLB = ISHLQ(4,I)

         IF (ISHLC.GT.0 .AND. ISHLD.GT.0 .AND.
     &       ISHLA.GT.0 .AND. ISHLB.GT.0) THEN

            IF (ISHLD .EQ. ISHLC) THEN
               NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
            ELSE
               NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
            END IF
            IF (ISHLB .EQ. ISHLA) THEN
               NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLA) + 1)/2
            ELSE
               NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
            END IF
            LINT1 = NUMCD*NUMAB ! actual space needed for (CD|AB)

C           Compute expected number of comparisons.
C           ---------------------------------------

            XPECTL = DBLE(LINT1)
            XPECT  = XPECT + XPECTL

C           Calculate shell quadruple (CD|AB).
C           ----------------------------------

            ISHLCD = CHO_F2SP(ITRI(ISHLC,ISHLD))
            ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
            IF (ISHLAB.LT.1 .OR. ISHLCD.LT.1) THEN
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF

            CALL CHO_DZERO(WORK(KINT1),LINT1)
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,WORK(KINT1),LINT1,.FALSE.)

C           Calculate integrals from Cholesky vectors.
C           ------------------------------------------

            CALL CHO_DBGINT_CHO(WORK(KINT1),NUMCD,NUMAB,WORK(KWRK),
     &                          LWRK/2,ERRMAX,ERRMIN,ERRRMS,NCMP,
     &                          ISHLCD,ISHLAB)

C           Write report.
C           -------------

            IF (NCMP .LT. 1) THEN
               WRITE(LUPRI,'(4(I5,1X),5X,A)')
     &         ISHLC,ISHLD,ISHLA,ISHLB,' !!! nothing compared !!! '
            ELSE
               XCMP  = DBLE(NCMP)
               XTCMP = XTCMP + XCMP
               RMS   = SQRT(ERRRMS/XCMP)
               IF (PRTLAB) THEN
                  CALL CHO_INTCHK_ID_OF(LABEL,I,-1)
                  WRITE(LUPRI,'(4(I5,1X),1P,3(D12.4,1X),A,A,A)')
     &            ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS,
     &            '(',LABEL,')'
               ELSE
                  WRITE(LUPRI,'(4(I5,1X),1P,3(D12.4,1X))')
     &            ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS
               END IF
            END IF

            IF (ABS(ERRMAX) .GT. ABS(GLMAX)) THEN
               GLMAX = ERRMAX
            END IF
            IF (ABS(ERRMIN) .LT. ABS(GLMIN)) THEN
               GLMIN = ERRMIN
            END IF
            GLRMS = GLRMS + ERRRMS

         END IF

      END DO

C     Print end of table.
C     -------------------

      WRITE(LUPRI,'(A)')
     & '--------------------------------------------------------------'
      IF (XTCMP .LT. 1.0D0) THEN
         WRITE(LUPRI,'(A,23X,A)')
     &   'Total:',' !!! nothing compared !!! '
      ELSE
         GLRMS = SQRT(GLRMS/XTCMP)
         WRITE(LUPRI,'(A,18X,1P,3(D12.4,1X))')
     &   'Total:',GLMIN,GLMAX,GLRMS
      END IF
      WRITE(LUPRI,'(A)')
     & '--------------------------------------------------------------'

C     Release all memory allocated here (and release seward memory).
C     --------------------------------------------------------------

      CALL XRLSMEM_INTS
      CALL GETMEM('DBGINT.2','FREE','REAL',KWRK,LWRK/2)
      CALL GETMEM('INTDBG.3','FLUSH','REAL',KINT1,LINTMX)
      CALL GETMEM('INTDBG.4','FREE','REAL',KINT1,LINTMX)

C     Check total number of comparisons.
C     ----------------------------------

      DO ISYM = 1,NSYM
         XLBAS(ISYM) = DBLE(NBAS(ISYM))
      END DO
      XNINT = 0.0D0
      DO ISYM = 1,NSYM
         XXLBST = 0.0D0
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .EQ. ISYMB) THEN
               XXLBST = XXLBST + XLBAS(ISYMA)*(XLBAS(ISYMA)+1.0D0)/2.0D0
            ELSE IF (ISYMA .GT. ISYMB) THEN
               XXLBST = XXLBST + XLBAS(ISYMA)*XLBAS(ISYMB)
            END IF
         END DO
         XNINT = XNINT + XXLBST*(XXLBST + 1.0D0)/2.0D0
      END DO

      IF (ABS(XTCMP-XPECT) .GT. 1.0D-15) THEN
         WRITE(LUPRI,'(/,A)')
     &   'WARNING: not all integrals checked:'
      ELSE
         WRITE(LUPRI,*)
      END IF
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number of integral comparisons    :',XTCMP
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number expected (full shell pairs):',XPECT
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number of unique integrals        :',XNINT

      END
