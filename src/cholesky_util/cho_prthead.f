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
      SUBROUTINE CHO_PRTHEAD(SKIPH)
C
C     Purpose: print Cholesky header.
C
      use ChoSubScr, only: Cho_SScreen, SSTau
#include "implicit.fh"
      LOGICAL SKIPH
#include "cholesky.fh"
#include "choprint.fh"

      PARAMETER (NADRMODE = 2, NALG = 6)

      CHARACTER*12 ALGORITHM(0:NALG)
      CHARACTER*13 ADRMODE(0:NADRMODE)
      CHARACTER*15 USED(2)

      DATA ALGORITHM /'     unknown',
     &                '    one-step','    two-step','       naive',
     &                'par one-step','par two-step','   par naive'/
      DATA ADRMODE /'      unknown','   word addr.','  direct acc.'/
      DATA USED /'(screening off)','(screening on) '/

      IF (LUPRI .LT. 1) THEN
         CALL CHO_QUIT('LUPRI undefined in Cholesky decomposition',101)
      END IF

      IF (.NOT.SKIPH) THEN
         WRITE(LUPRI,'(//,80A)') ('*', I=1,80)
         WRITE(LUPRI,'(A,78X,A)') ('*',I=1,2)
         WRITE(LUPRI,'(A,10X,A,10X,A)')
     & '*','Cholesky Decomposition of Two-Electron Repulsion Integrals',
     &   '*'
         WRITE(LUPRI,'(A,78X,A)') ('*',I=1,2)
         WRITE(LUPRI,'(80A)') ('*', I=1,80)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)

         IF (RSTDIA) THEN
            WRITE(LUPRI,'(/,A)')
     &      '***** Using Restart Integral Diagonal *****'
         END IF
         IF (RSTCHO) THEN
            IF (RSTDIA) THEN
               WRITE(LUPRI,'(A)')
     &         '***** Using Restart Cholesky Vectors  *****'
            ELSE
               WRITE(LUPRI,'(/,A)')
     &         '***** Using Restart Cholesky Vectors  *****'
            END IF
         END IF
      END IF

      IF (IPRINT .GE. INF_INIT) THEN
         IF (SCDIAG) THEN
            IUSE = 2
         ELSE
            IUSE = 1
         END IF
         IF (CHO_DECALG.LT.1 .OR. CHO_DECALG.GT.NALG) THEN
            IALG = 0
         ELSE
            IALG = CHO_DECALG
         END IF
         IF (.NOT. SKIPH) CALL CHO_HEAD('Configuration','=',80,LUPRI)
         WRITE(LUPRI,'(A,A)')
     &   'Decomposition algorithm                   : ',ALGORITHM(IALG)
         IF (CHO_1CENTER) THEN
            WRITE(LUPRI,'(A)')
     &      '1-center decomposition                    :          Yes'
            IF (CHO_NO2CENTER) THEN
               WRITE(LUPRI,'(A)')
     &        'Exclusion of 2-center diagonals           :          Yes'
            ELSE
               WRITE(LUPRI,'(A)')
     &        'Exclusion of 2-center diagonals           :           No'
            END IF
            IF (CHO_SIMRI) THEN
               WRITE(LUPRI,'(A,1P,D12.4)')
     &         'Simulation of RI, threshold               : ',THR_SIMRI
            END IF
         ELSE
            WRITE(LUPRI,'(A)')
     &      '1-center decomposition                    :           No'
         END IF
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Decomposition threshold                   : ',THRCOM
         IF (CHO_PRESCREEN) THEN
            WRITE(LUPRI,'(A,1P,D12.4)')
     &      'Initial diagonal prescreening             : ',THR_PRESCREEN
         END IF
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Initial diagonal screening                : ',THRDIAG
         WRITE(LUPRI,'(A,1P,D12.4,1X,A)')
     &   'First  screening damping                  : ',DAMP(1),
     &                                                  USED(IUSE)
         WRITE(LUPRI,'(A,1P,D12.4,1X,A)')
     &   'Second screening damping                  : ',DAMP(2),
     &                                                  USED(IUSE)
         IF (CHO_USEABS) THEN
            WRITE(LUPRI,'(A)')
     &      'Absolute values used in diagonal screening:          Yes'
         ELSE
            WRITE(LUPRI,'(A)')
     &      'Absolute values used in diagonal screening:           No'
         END IF
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Threshold for negative  diagonal zeroing  : ',THRNEG
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Threshold for warning about neg. diagonal : ',WARNEG
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Threshold for too negative diagonal       : ',TOONEG
         WRITE(LUPRI,'(A,1P,D12.4)')
     &   'Span factor                               : ',SPAN
         WRITE(LUPRI,'(A,I12)')
     &   'Max. #Cholesky vectors per symmetry       : ',MAXVEC
         WRITE(LUPRI,'(A,I12)')
     &   'Max. #reduced sets (i.e., integral passes): ',MAXRED
         WRITE(LUPRI,'(A,I12)')
     &   'Min. #qualified required for decomposition: ',MINQUAL
         WRITE(LUPRI,'(A,I12)')
     &   'Max. #qualified per symmetry              : ',MAXQUAL
         IF (N2_QUAL .EQ. 0) THEN
            XF = -9.99999999D15
         ELSE
            X1 = DBLE(N1_QUAL)
            X2 = DBLE(N2_QUAL)
            XF = 1.0D2*X1/X2
         END IF
         WRITE(LUPRI,'(A,5X,F7.4,A)')
     &   'Max. memory fraction for qualified columns: ',XF,'%'
         IF (MXSHPR .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     &      'Max. #shell pair allowed per integral pass:      generic'
         ELSE
            WRITE(LUPRI,'(A,I12)')
     &      'Max. #shell pair allowed per integral pass: ',MXSHPR
         END IF
         IF (IALQUA .EQ. 0) THEN
            WRITE(LUPRI,'(A)')
     &      'Qualification algorithm                   : dalton-style'
         ELSE IF (IALQUA .EQ. 1) THEN
            WRITE(LUPRI,'(A)')
     &      'Qualification algorithm                   :   sequential'
         ELSE
            WRITE(LUPRI,'(A)')
     &      'Qualification algorithm                   :      sorting'
         END IF
         IF (CHO_IOVEC .EQ. 1) THEN
            WRITE(LUPRI,'(A)')
     &      'Algorithm for Cholesky vector I/O         :  rs2rs/batch'
         ELSE IF (CHO_IOVEC .EQ. 2) THEN
            WRITE(LUPRI,'(A)')
     &      'Algorithm for Cholesky vector I/O         : buffer/rs2rs'
         ELSE IF (CHO_IOVEC .EQ. 3) THEN
            WRITE(LUPRI,'(A)')
     &      'Algorithm for Cholesky vector I/O         : lrgbuf/rs2rs'
         ELSE IF (CHO_IOVEC .EQ. 4) THEN
            WRITE(LUPRI,'(A)')
     &      'Algorithm for Cholesky vector I/O         : fxdbuf/rs2rs'
         ELSE
            WRITE(LUPRI,'(A)')
     &      'Algorithm for Cholesky vector I/O         : copy via rs1'
         END IF
         IADRMODE = MAX(MIN(CHO_ADRVEC,NADRMODE),0)
         WRITE(LUPRI,'(A,A13)')
     &   'Address mode for Cholesky vector I/O      : ',
     &   ADRMODE(IADRMODE)
         XF2 = 1.0D2*FRAC_CHVBUF
         WRITE(LUPRI,'(A,5X,F7.4,A)')
     &   'Memory fraction used as vector buffer     : ',XF2,'%'
         IF (CHO_SSCREEN) THEN
            WRITE(LUPRI,'(A,1P,D12.4)')
     &      'Screening threshold for vector subtraction: ',SSTAU
         END IF
         IF (CHO_DECALG .EQ. 5) THEN
            WRITE(LUPRI,'(A,I12)')
     &      'Block size (blocked Z vector array)       : ',BLOCKSIZE
         END IF
      END IF

      END
