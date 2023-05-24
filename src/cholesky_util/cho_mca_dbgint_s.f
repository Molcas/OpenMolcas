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
      SUBROUTINE CHO_MCA_DBGINT_S(ISHLQ,NSHLQ,PRTLAB)
!
!     Purpose: regenerate and check a specified set of integrals.
!
!     NOTE: this is *not* meant for production calculations, only for
!           debugging, as:
!           1) full Cholesky vectors are read
!           2) calculations are performed in full (no use of red. sets
!              apart from first)
!
      use ChoArr, only: nBstSh
      Implicit Real*8 (a-h,o-z)
      INTEGER ISHLQ(4,NSHLQ)
      LOGICAL PRTLAB
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      CHARACTER(LEN=16), PARAMETER:: SECNAM = 'CHO_MCA_DBGINT_S'

      INTEGER, EXTERNAL:: CHO_F2SP

      Real*8 XLBAS(8)

      CHARACTER(LEN=8) LABEL

      Real*8, Allocatable:: INT1(:), WRK(:)

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      ITRI(I,J)=MAX(I,J)*(MAX(I,J)-3)/2+I+J

!     Return if nothing specified.
!     ----------------------------

      IF (NSHLQ .LT. 1) RETURN

!     Force computation of full shell quadruple.
!     ------------------------------------------

      IF (IFCSEW .NE. 1) THEN
         WRITE(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',
     &                  IFCSEW,' to 1.'
         WRITE(LUPRI,*) SECNAM,
     &   ': memory demands are significantly increased by this!'
         IFCSEW = 1
      END IF

!     Initializations.
!     ----------------

      GLMAX = 0.0D0
      GLMIN = 1.0D15
      GLRMS = 0.0D0
      XTCMP = 0.0D0
      XPECT = 0.0D0

!     Make first reduced set the current reduced set.
!     -----------------------------------------------

      CALL CHO_RSCOPY(1,2)

!     Allocate memory for largest integral quadruple.
!     -----------------------------------------------

      LINTMX = MX2SH*MX2SH
      Call mma_allocate(INT1,LINTMX,Label='INT1')

!     Allocate max. memory
!     ----------------------------------------------------------

      Call mma_maxDBLE(LWRK)
      Call mma_allocate(WRK,LWRK/2,Label='WRK')
      CALL XSETMEM_INTS(LWRK/2)

!     Print header.
!     -------------

      CALL CHO_HEAD('Integral Error Analysis','=',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A)')
     & '    C     D     A     B   Abs. Min.    Abs. Max.      RMS',
     & '--------------------------------------------------------------'

!     Loop over specified shell quadruples.
!     -------------------------------------

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

!           Compute expected number of comparisons.
!           ---------------------------------------

            XPECTL = DBLE(LINT1)
            XPECT  = XPECT + XPECTL

!           Calculate shell quadruple (CD|AB).
!           ----------------------------------

            ISHLCD = CHO_F2SP(ITRI(ISHLC,ISHLD))
            ISHLAB = CHO_F2SP(ITRI(ISHLA,ISHLB))
            IF (ISHLAB.LT.1 .OR. ISHLCD.LT.1) THEN
               CALL CHO_QUIT('CHO_F2SP<1 in '//SECNAM,103)
            END IF

            INT1(1:LINT1)=0.0d0
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,INT1,LINT1,.FALSE.)

!           Calculate integrals from Cholesky vectors.
!           ------------------------------------------

            CALL CHO_DBGINT_CHO(INT1,NUMCD,NUMAB,WRK,
     &                          LWRK/2,ERRMAX,ERRMIN,ERRRMS,NCMP,
     &                          ISHLCD,ISHLAB)

!           Write report.
!           -------------

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

!     Print end of table.
!     -------------------

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

!     Release all memory allocated here (and release seward memory).
!     --------------------------------------------------------------

      CALL XRLSMEM_INTS
      Call mma_deallocate(WRK)
      Call mma_deallocate(INT1)

!     Check total number of comparisons.
!     ----------------------------------

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
