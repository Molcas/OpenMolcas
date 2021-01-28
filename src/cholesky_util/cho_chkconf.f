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
      SUBROUTINE CHO_CHKCONF(NCONFL,VERBOSE)
C
C     Purpose: check configuration, return the number of errors NCONFL.
C
      use ChoSubScr, only: Cho_SScreen, SSTau
      IMPLICIT NONE
      INTEGER NCONFL
      LOGICAL VERBOSE
#include "cholesky.fh"
#include "choorb.fh"
#include "chosimri.fh"
#include "cho_para_info.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_CHKCONF')

      LOGICAL REPORT
      INTEGER NNN, MMM, INEGRR
      REAL*8  XLBUF, XMBUF

C     Initialize.
C     -----------

      NCONFL = 0
      REPORT = VERBOSE

C     Check that output unit is appropriately set.
C     (Upper bound not checked, as it may vary.)
C     --------------------------------------------

      IF (LUPRI .LT. 1) THEN
         NCONFL = NCONFL + 1
         REPORT = .FALSE.
      END IF

C     Check decomposition algorithm.
C     ------------------------------

      IF (CHO_DECALG.LT.1 .OR. CHO_DECALG.GT.CHO_NDECALG) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,I4)') 'Illegal decomposition algorithm, ',
     &                     'CHO_DECALG = ',CHO_DECALG
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (Cho_Real_Par) THEN
         IF (CHO_DECALG.NE.4 .AND. CHO_DECALG.NE.5 .AND.
     &       CHO_DECALG.NE.6) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,I4)')
     &                        'Illegal decomposition algorithm, ',
     &                        'CHO_DECALG = ',CHO_DECALG
               WRITE(LUPRI,'(A,A)')
     &                        'Only parallel algorithm is allowed ',
     &                        'for parallel execution'
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Cancel exclusion of 2-center diagonals in case of symmetry.
C     ----------------------------------------------------------

      IF (CHO_NO2CENTER) THEN
         IF (NSYM .NE. 1) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') 'Exclusion of 2-center diagonals ',
     &                        'only implemented for C1 point group.'
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Cancel 1-center decomposition in case of symmetry.
C     --------------------------------------------------

      IF (CHO_1CENTER) THEN
         IF (NSYM .NE. 1) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') '1-center decomposition ',
     &                        'only implemented for C1 point group.'
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Checks specific to RI simulation.
C     ---------------------------------

      IF (CHO_SIMRI) THEN
         IF (.NOT. CHO_1CENTER) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') '1-center decomposition required ',
     &                        'for RI simulation.'
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (CHO_DECALG .NE. 2) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') 'RI simulation can only be ',
     &                        'executed with the two-step algorithm.'
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (RSTCHO) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') 'RI simulation can not be ',
     &                        'executed with decomposition restart.'
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Check for conflicts for specific algorithms.
C     --------------------------------------------

      IF (CHO_DECALG.EQ.2 .OR. CHO_DECALG.EQ.5) THEN
         IF (.NOT. SCDIAG) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)')
     &                        'Screening must be turned on for two-',
     &                        'step decomposition algorithm.'
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (IFCSEW .NE. 2) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,A)') 'The interface to Seward must ',
     &                        'be "2" (reduced set communication)',
     &                        ' for two-step algorithm.'
               WRITE(LUPRI,'(A,I4)') 'Current value is IFCSEW=',IFCSEW
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (RSTCHO) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A)') 'Decomposition restart is not ',
     &                              'possible for two-step algorithm.'
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

      IF (CHO_DECALG.EQ.5) THEN
         IF (BLOCKSIZE.LT.1) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,I8)')
     &         'BlockSize must be strictly positive.',
     &         'Current value is BlockSize=',BlockSize
            END IF
            NCONFL=NCONFL+1
         END IF
      END IF

C     Check max. number of Cholesky vectors and reduced sets.
C     -------------------------------------------------------

      IF (MAXVEC .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,I8)') 'Max. number of vectors < 1: ',MAXVEC
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (MAXRED .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,I8)')
     &      'Max. number of reduced sets < 1: ',MAXRED
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Check qualification algorithm.
C     ------------------------------

      IF (IALQUA .LT. 0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,I4,A)')
     &                     'Qualification algorithm reset from ',
     &                     IALQUA,' to 1'
         END IF
         IALQUA = 1
      ELSE IF (IALQUA .GT. NALQUA) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &                   'Qualification algorithm reset from ',
     &                   IALQUA,' to ',NALQUA
         END IF
         IALQUA = NALQUA
      END IF

C     Decomposition threshold.
C     ------------------------

      IF (THRCOM .LT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6,A,D15.6,A)')
     &                     'Decomposition threshold not positive: ',
     &                     THRCOM,' (default value: ',THRDEF,
     &                     ')'
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Diagonal prescreening threshold.
C     --------------------------------

      IF (CHO_PRESCREEN) THEN
         IF (THR_PRESCREEN .GT. THRCOM) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,1P,D15.6,A,D15.6)')
     &                        'Diagonal prescreening threshold is ',
     &                        'greater than decomposition threshold: ',
     &                        THR_PRESCREEN,' > ',THRCOM
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Screening threshold for vector subtraction.
C     -------------------------------------------

      IF (CHO_SSCREEN) THEN
         IF (SSTAU .LT. 0.0D0) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,1P,D15.6)')
     &                        'Screening threshold for vector ',
     &                        'subtraction not positive: ',SSTAU
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Buffer length.
C     --------------

      IF (LBUF .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,I8)') 'Buffer length < 1: ',LBUF
         END IF
         NCONFL = NCONFL + 1
      ELSE
         XLBUF = DBLE(LBUF)
         XMBUF = DBLE(NBAST)*(DBLE(NBAST)+1.0D0)/2.0D0
         IF (XLBUF .GT. XMBUF) THEN ! make sure LBUF is not too large
            LBUF = NINT(XMBUF)
         END IF
      END IF

C     Memory split for qualified columns.
C     -----------------------------------

      IF (N1_QUAL .GE. N2_QUAL) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,2I5)') 'N1_QUAL >= N2_QUAL: ',
     &                     N1_QUAL,N2_QUAL
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (N1_QUAL.LT.1 .OR. N2_QUAL.LT.1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,2I5)')
     &                     'N1_QUAL and/or N2_QUAL < 1: ',
     &                     N1_QUAL,N2_QUAL
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Memory split for buffered reading of previous vectors.
C     Max. #vectors in each call to DGEMM in subtraction.
C     (CHO_IOVEC=3,4 only.)
C     ------------------------------------------------------

      IF (CHO_IOVEC.EQ.3 .OR. CHO_IOVEC.EQ.4) THEN
         IF (N2_VECRD .LT. N1_VECRD) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,2I5)') 'N1_VECRD >= N2_VECRD: ',
     &                        N1_VECRD,N2_VECRD
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (N1_VECRD.LT.1 .OR. N2_VECRD.LT.1) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,2I5)') 'N1_VECRD and/or N2_VECRD < 1: ',
     &                        N1_VECRD,N2_VECRD
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (N_SUBTR .LT. 1) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,I8)') 'N_SUBTR: ',N_SUBTR
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

C     Memory fraction used for vector buffer.
C     ---------------------------------------

      IF (FRAC_CHVBUF .LT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6,A)') 'FRAC_CHVBUF=',FRAC_CHVBUF,
     &                     ' resetting value to 0.0D0'
         END IF
         FRAC_CHVBUF = 0.0D0
      ELSE IF (FRAC_CHVBUF .GT. 0.9D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6,A)') 'FRAC_CHVBUF=',FRAC_CHVBUF,
     &                     ' resetting value to 0.9D0'
         END IF
         FRAC_CHVBUF = 0.9D0
      END IF

C     Threshold for discarding elements of initial diagonal.
C     ------------------------------------------------------

      IF (THRDIAG .GT. THRCOM) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,1P,D15.6)')
     &                     'Threshold for discarding initial diagonals',
     &                     ': ',THRDIAG
            WRITE(LUPRI,'(A,1P,D15.6)')
     &                     'is larger than decomposition threshold: ',
     &                     THRCOM
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Damping factors.
C     ----------------

      IF (DAMP(1) .LT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6)')
     &      'First damping factor < 1: ',DAMP(1)
         END IF
         NCONFL = NCONFL + 1
      END IF
      IF (DAMP(2) .LT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6)')
     &      'Second damping factor < 1: ',DAMP(2)
         END IF
         NCONFL = NCONFL + 1
      END IF

C     Span factor.
C     ------------

      IF (SPAN .GT. 1.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6)') 'Span factor > 1: ',SPAN
         END IF
         NCONFL = NCONFL + 1
      ELSE IF (ABS(1.0D0-SPAN).lt.1.0d-4) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A)')
     &      'Span factor is too close to 1. Will use 0.9999 instead.'
         END IF
         SPAN=0.9999d0
      END IF

C     Max. #shell pairs allowed before proceeding to deco.
C     Max. #qualifieds allowed to proceed to decomposition.
C     Min. #qualifieds needed to proceed to decomposition.
C     -----------------------------------------------------

      NNN = NSHELL*(NSHELL + 1)/2
      IF (MXSHPR.LT.0 .OR. MXSHPR.GT.NNN) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,I8)') 'Max. #shell pairs allowed ',
     &                     'before proceeding to deco. is: ',MXSHPR
            WRITE(LUPRI,'(A)') 'Resetting generic: 0'
         END IF
         MXSHPR = 0
      END IF
      IF (MAXQUAL .LT. 1) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,I8)') 'Max. number of qualifieds is ',
     &                     'non-positive: ',MAXQUAL
            WRITE(LUPRI,'(A,A,I8)') 'Using instead: ',ABS(MAXQUAL)
         END IF
         MAXQUAL = ABS(MAXQUAL)
      END IF
      MMM = NSYM*MAXQUAL
      IF (MINQUAL.LT.1 .OR. MINQUAL.GT.MMM) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,I8)') 'Min. #qualified needed ',
     &                     'before proceeding to deco. is: ',MINQUAL
            WRITE(LUPRI,'(A,I8)') 'resetting to: ',MMM
         END IF
         MINQUAL = MMM
      END IF

C     Handling of negative diagonals.
C     -------------------------------

      INEGRR = 0
      IF (THRNEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,1P,D15.6)')
     &                     'Threshold for zeroing neg. diag. > 0: ',
     &                     THRNEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (WARNEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,1P,D15.6)')
     &                     'Threshold for warning about neg. diag. ',
     &                     ' > 0: ',WARNEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (TOONEG .GT. 0.0D0) THEN
         IF (REPORT) THEN
            WRITE(LUPRI,'(A,A,1P,D15.6)')
     &                      'Threshold for shutdown due to neg. diag. ',
     &                     ' > 0: ',TOONEG
         END IF
         INEGRR = INEGRR + 1
         NCONFL = NCONFL + 1
      END IF
      IF (INEGRR .EQ. 0) THEN
         IF (THRNEG .LT. WARNEG) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,1P,D15.6,D15.6)')
     &                      'Threshold for zeroing neg. diag. > ',
     &                       'threshold for warning about neg. diag.: ',
     &                        THRNEG,WARNEG
            END IF
            NCONFL = NCONFL + 1
         END IF
         IF (WARNEG .LT. TOONEG) THEN
            IF (REPORT) THEN
               WRITE(LUPRI,'(A,A,1P,D15.6,D15.6)')
     &                     'Threshold for warning about neg. diag. > ',
     &                     'threshold for shutdown due to neg. diag.: ',
     &                     WARNEG,TOONEG
            END IF
            NCONFL = NCONFL + 1
         END IF
      END IF

      END
