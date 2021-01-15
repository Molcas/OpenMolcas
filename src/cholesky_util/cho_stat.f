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
      SUBROUTINE CHO_STAT()
C
C     Purpose: print statistics from decomposition.
C
      USE Para_Info, ONLY: nProcs, Is_Real_Par
      use ChoArr, only: nDimRS, IntMap
#include "implicit.fh"
#include "cholesky.fh"
#include "choprint.fh"
#include "choorb.fh"
#include "choptr.fh"
#include "chosubscr.fh"
#include "WrkSpc.fh"

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_STAT')

      CHARACTER*25 STRING
      CHARACTER*2  UNT

      PARAMETER (NTAU = 5, DTAU = 1.0D-1)
      REAL*8 TAU(NTAU), XC(NTAU)

      REAL*8 XXBST(8), VCSTOR(8)

      LOGICAL DOCPCT, DOWPCT, CHO_SSCREEN_SAVE, PARALG

      PARAMETER (N2 = INFVEC_N2)

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      NNBSTRSH(I,J,K)=IWORK(ip_NNBSTRSH-1+NSYM*NNSHL*(K-1)+NSYM*(J-1)+I)
      DSPNM(I)=WORK(ip_DSPNM-1+I)


      PARALG=CHO_DECALG.EQ.4 .OR. CHO_DECALG.EQ.5 .OR. CHO_DECALG.EQ.6

C     Overall header.
C     ---------------

      CALL CHO_HEAD('Cholesky Decomposition Statistics','=',80,LUPRI)
      IF (RSTDIA) THEN
         WRITE(LUPRI,'(/,A)')
     &   'Calculation restarted from diagonal on disk'
      END IF
      IF (RSTCHO) THEN
         IF (RSTDIA) THEN
            WRITE(LUPRI,'(A)')
     &      'Calculation restarted from Cholesky vectors on disk'
         ELSE
            WRITE(LUPRI,'(/,A)')
     &      'Calculation restarted from Cholesky vectors on disk'
         END IF
      END IF

C     Configuration.
C     --------------

      CALL CHO_HEAD('Configuration','-',80,LUPRI)
      IPRSAV = IPRINT
      IPRINT = INF_INIT + 1
      CALL CHO_PRTHEAD(.TRUE.)
      IPRINT = IPRSAV

C     Vector statistics.
C     ------------------

      CALL CHO_HEAD('Vector statistics','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,A,/,A,A)')
     & 'Sym.        N      Full      Mmax         M    M/Full',
     & '    M/Mmax       M/N',
     & '-----------------------------------------------------',
     & '---------------------'
      NTOT = 0
      DO ISYM = 1,NSYM
         NTOT = NTOT + NUMCHO(ISYM)
         NN = 0
         XXBST(ISYM) = 0.0D0
         DO JSYM = 1,NSYM
            KSYM = MULD2H(JSYM,ISYM)
            IF (JSYM .GT. KSYM) THEN
               NN = NN + NBAS(JSYM)*NBAS(KSYM)
               XXBST(ISYM) = XXBST(ISYM)
     &         + DBLE(NBAS(JSYM))*DBLE(NBAS(KSYM))
            ELSE IF (JSYM .EQ. KSYM) THEN
               NN = NN + NBAS(JSYM)*(NBAS(JSYM) + 1)/2
               XXBST(ISYM) = XXBST(ISYM)
     &         + DBLE(NBAS(JSYM))*(DBLE(NBAS(JSYM)) + 1.0D0)/2.0D0
            END IF
         END DO
         IF (NN .GT. 0) THEN
            YYY = DBLE(NUMCHO(ISYM))/XXBST(ISYM)
         ELSE
            YYY = 9.0D9
         END IF
         IF (NNBSTR(ISYM,1) .GT. 0) THEN
            XX = DBLE(NNBSTR(ISYM,1))
            YY = DBLE(NUMCHO(ISYM))/XX
         ELSE
            YY = 9.0D9
         END IF
         IF (NBAS(ISYM) .NE. 0) THEN
            X = DBLE(NBAS(ISYM))
            Y = DBLE(NUMCHO(ISYM))/X
         ELSE
            Y = 9.0D9
         END IF
         WRITE(LUPRI,'(I3,4(1X,I9),3(1X,F9.4))')
     &   ISYM,NBAS(ISYM),NN,NNBSTR(ISYM,1),NUMCHO(ISYM),YYY,YY,Y
      END DO
      WRITE(LUPRI,'(A,A)')
     & '-----------------------------------------------------',
     & '---------------------'
      NN = NBAST*(NBAST + 1)/2
      IF (NN .GT. 0) THEN
         XXX = DBLE(NN)
         YYY = DBLE(NUMCHT)/XXX
      ELSE
         YYY = 9.0D9
      END IF
      IF (NNBSTRT(1) .GT. 0) THEN
         XX = DBLE(NNBSTRT(1))
         YY = DBLE(NUMCHT)/XX
      ELSE
         YY = 9.0D9
      END IF
      IF (NBAST .NE. 0) THEN
         X = DBLE(NBAST)
         Y = DBLE(NUMCHT)/X
      ELSE
         Y = 9.0D9
      END IF
      WRITE(LUPRI,'(3X,4(1X,I9),3(1X,F9.4))')
     & NBAST,NN,NNBSTRT(1),NUMCHT,YYY,YY,Y
      WRITE(LUPRI,'(A,A)')
     & '-----------------------------------------------------',
     & '---------------------'
      IF (NTOT .NE. NUMCHT) THEN
         WRITE(LUPRI,'(A)')
     &   'WARNING: total number of vectors is wrong!!!'
      END IF

      CALL CHO_GETSTOR(VCSTOR)

      WRITE(LUPRI,'(/,A,/,A,/,A)')
     & '                       %Saving relative to',
     & 'Sym.     Storage      1st Red. Set     Full',
     & '---------------------------------------------'
      VCTOT = 0.0D0
      X1TOT = 0.0D0
      X2TOT = 0.0D0
      XXTOT = 0.0D0
      DO ISYM = 1,NSYM
         X1DIM = DBLE(NNBSTR(ISYM,1))
         X1RED = X1DIM*DBLE(NUMCHO(ISYM))
         IF (X1RED .GT. 0.0D0) THEN
            SAV1  = 1.0D2*(X1RED - VCSTOR(ISYM))/X1RED
         ELSE
            SAV1 = 9.0D9
         END IF
         XX = XXBST(ISYM)*DBLE(NUMCHO(ISYM))
         IF (XX .GT. 0.0D0) THEN
            SAV2  = 1.0D2*(XX - VCSTOR(ISYM))/XX
         ELSE
            SAV2 = 9.0D9
         END IF
         XGB = VCSTOR(ISYM)*8.0D0/1.024D3
         UNT = 'kb'
         IF (XGB .GT. 1.0D3) THEN
            XGB = XGB/1.024D3
            UNT = 'Mb'
            IF (XGB .GT. 1.0D3) THEN
               XGB = XGB/1.024D3
               UNT = 'Gb'
               IF (XGB .GT. 1.0D3) THEN
                  XGB = XGB/1.024D3
                  UNT = 'Tb'
               END IF
            END IF
         END IF
         WRITE(LUPRI,'(I2,4X,F10.3,1X,A,4X,F9.4,4X,F9.4)')
     &   ISYM,XGB,UNT,SAV1,SAV2
         VCTOT = VCTOT + VCSTOR(ISYM)
         X1TOT = X1TOT + X1RED
         X2TOT = X2TOT + X1DIM*(X1DIM + 1.0D0)/2.0D0
         XXTOT = XXTOT + XX
      END DO
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------'
      IF (X1TOT .GT. 0.0D0) THEN
         SAV1  = 1.0D2*(X1TOT - VCTOT)/X1TOT
      ELSE
         SAV1 = 9.0D9
      END IF
      IF (XXTOT .GT. 0.0D0) THEN
         SAV2  = 1.0D2*(XXTOT - VCTOT)/XXTOT
      ELSE
         SAV2 = 9.0D9
      END IF
      XGB = VCTOT*8.0D0/1.024D3
      UNT = 'kb'
      IF (XGB .GT. 1.0D3) THEN
         XGB = XGB/1.024D3
         UNT = 'Mb'
         IF (XGB .GT. 1.0D3) THEN
            XGB = XGB/1.024D3
            UNT = 'Gb'
            IF (XGB .GT. 1.0D3) THEN
               XGB = XGB/1.024D3
               UNT = 'Tb'
            END IF
         END IF
      END IF
      WRITE(LUPRI,'(A6,F10.3,1X,A,4X,F9.4,4X,F9.4)')
     & 'Total:',XGB,UNT,SAV1,SAV2
      WRITE(LUPRI,'(A)')
     & '---------------------------------------------'
      PCT = 1.0D2*VCTOT/X2TOT
      WRITE(LUPRI,'(A,F11.6,A)')
     & 'Total storage corresponds to',PCT,
     & '% of the 1st reduced set integral matrix.'

      CALL CHO_STAT_PARENTDIAG()

C     Integral statistics.
C     --------------------

      MAXCAL = 0
      NCAL   = 0
      NREP   = 0
      NCALL  = 0
      DO ISHLAB = 1,NNSHL
         IF (INTMAP(ISHLAB) .GT. 0) THEN
            MAXCAL = MAX(MAXCAL,INTMAP(ISHLAB))
            NCAL   = NCAL  + 1
            NCALL  = NCALL + INTMAP(ISHLAB)
            IF (INTMAP(ISHLAB) .GT. 1) NREP = NREP + 1
         END IF
      END DO

      XXSHL = DBLE(NNSHL)
      XCAL  = DBLE(NCAL)
      XREP  = DBLE(NREP)

      CALL CHO_HEAD('Integral statistics','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,I10)')
     & '#Shells                  :',NSHELL
      WRITE(LUPRI,'(A,I10)')
     & '#Shell Pair Distributions:',NNSHL
      WRITE(LUPRI,'(A,I10)')
     & '#Integral passes         :',XNPASS
      WRITE(LUPRI,'(A,I10)')
     & '#Calls to integral prog. :',NCALL
      WRITE(LUPRI,'(A,I10,A,F8.3,A,A)')
     & '#Shell Pairs Calculated  :',NCAL,' (',XCAL*1.0D2/XXSHL,' %',
     & ' of total)'
      WRITE(LUPRI,'(A,I10,A,F8.3,A,A)')
     & '#Shell Pairs Repeated    :',NREP,' (',XREP*1.0D2/XCAL,' %',
     & ' of calculated)'

      WRITE(LUPRI,'(/,A,/,A)')
     & '#Calculations     #Shell Pairs   Percentage',
     & '-------------------------------------------'
      DO ICAL = 1,MAXCAL
         N = 0
         DO ISHLAB = 1,NNSHL
            IF (INTMAP(ISHLAB) .EQ. ICAL) N = N + 1
         END DO
         IF (N .GT. 0) THEN
            X = DBLE(N)
            WRITE(LUPRI,'(I12,6X,I12,5X,F8.3)') ICAL,N,X*1.0D2/XXSHL
         END IF
      END DO
      WRITE(LUPRI,'(A)')
     & '-------------------------------------------'

C     Section timings.
C     ----------------

      IF (IPRINT .GE. INF_TIMING) THEN

      TCINI = TIMSEC(2,1) - TIMSEC(1,1)
      TWINI = TIMSEC(4,1) - TIMSEC(3,1)
      TCDIA = TIMSEC(2,2) - TIMSEC(1,2)
      TWDIA = TIMSEC(4,2) - TIMSEC(3,2)
      TCDEC = TIMSEC(2,3) - TIMSEC(1,3)
      TWDEC = TIMSEC(4,3) - TIMSEC(3,3)
      TCCHD = TIMSEC(2,4) - TIMSEC(1,4)
      TWCHD = TIMSEC(4,4) - TIMSEC(3,4)
      TCCHA = TIMSEC(2,5) - TIMSEC(1,5)
      TWCHA = TIMSEC(4,5) - TIMSEC(3,5)
      TCREO = TIMSEC(2,6) - TIMSEC(1,6)
      TWREO = TIMSEC(4,6) - TIMSEC(3,6)
      TCDIS = TIMSEC(2,7) - TIMSEC(1,7)
      TWDIS = TIMSEC(4,7) - TIMSEC(3,7)
      TCFIN = TIMSEC(2,8) - TIMSEC(1,8)
      TWFIN = TIMSEC(4,8) - TIMSEC(3,8)

      CALL CHO_HEAD('Section timings','-',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A,/,A)')
     &'                                 CPU time          Wall time',
     &'Section                      hours min. sec.    hours min. sec.',
     &'---------------------------------------------------------------'
      STRING = 'Initialization           '
      CALL CHO_CNVTIM(TCINI,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWINI,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Diagonal setup           '
      CALL CHO_CNVTIM(TCDIA,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWDIA,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Cholesky decomposition   '
      CALL CHO_CNVTIM(TCDEC,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWDEC,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      STRING = 'Diagonal check           '
      CALL CHO_CNVTIM(TCCHD,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWCHD,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      IF (CHO_INTCHK) THEN
         STRING = 'Integral check (debug)   '
         CALL CHO_CNVTIM(TCCHA,IHC,IMC,SCC)
         CALL CHO_CNVTIM(TWCHA,IHW,IMW,SCW)
         WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     &   STRING,IHC,IMC,SCC,IHW,IMW,SCW
      END IF
      IF (CHO_REORD) THEN
         STRING = 'Vector reordering        '
         CALL CHO_CNVTIM(TCREO,IHC,IMC,SCC)
         CALL CHO_CNVTIM(TWREO,IHW,IMW,SCW)
         WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     &    STRING,IHC,IMC,SCC,IHW,IMW,SCW
      END IF
      IF (CHO_FAKE_PAR .AND. NPROCS.GT.1 .AND. Is_Real_Par()) THEN
         STRING = 'Vector distribution      '
         CALL CHO_CNVTIM(TCDIS,IHC,IMC,SCC)
         CALL CHO_CNVTIM(TWDIS,IHW,IMW,SCW)
         WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     &    STRING,IHC,IMC,SCC,IHW,IMW,SCW
      END IF
      STRING = 'Finalization             '
      CALL CHO_CNVTIM(TCFIN,IHC,IMC,SCC)
      CALL CHO_CNVTIM(TWFIN,IHW,IMW,SCW)
      WRITE(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)')
     & STRING,IHC,IMC,SCC,IHW,IMW,SCW
      WRITE(LUPRI,'(A)')
     &'---------------------------------------------------------------'

      END IF ! IPRINT .GT. INF_TIMING

C     Timing of decomposition driver.
C     -------------------------------

      IF (DID_DECDRV) THEN

         CMISC = TDECDRV(1)
         WMISC = TDECDRV(2)
         DO J = 1,NINTEG
            CMISC = CMISC - TINTEG(1,J)
            WMISC = WMISC - TINTEG(2,J)
         END DO
         DO J = 1,NDECOM
            CMISC = CMISC - TDECOM(1,J)
            WMISC = WMISC - TDECOM(2,J)
         END DO
         DO J = 1,NMISC
            CMISC = CMISC - TMISC(1,J)
            WMISC = WMISC - TMISC(2,J)
         END DO
         CPCT = -9.0D9
         CFAC = -9.0D9
         WPCT = -9.0D9
         WFAC = -9.0D9
         IF (TDECDRV(1) .GT. 0.0D0) THEN
            DOCPCT = .TRUE.
            CFAC   = 1.0D2/TDECDRV(1)
         ELSE
            DOCPCT = .FALSE.
         END IF
         IF (TDECDRV(2) .GT. 0.0D0) THEN
            DOWPCT = .TRUE.
            WFAC   = 1.0D2/TDECDRV(2)
         ELSE
            DOWPCT = .FALSE.
         END IF

         CALL CHO_HEAD('Timing of decomposition driver','-',80,LUPRI)
         WRITE(LUPRI,'(/,A,A,/,A,A)')
     &    'Task           Component           CPU (min.)     %',
     &    '   Wall (min.)     %',
     &    '---------------------------------------------------',
     &    '----------------------'
         CTIM = TINTEG(1,1)
         WTIM = TINTEG(2,1)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Integrals      calculation        ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TINTEG(1,2)
         WTIM = TINTEG(2,2)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               I/O, qualifieds    ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,1)
         WTIM = TDECOM(2,1)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Decomposition  I/O, qualifieds    ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,2)
         WTIM = TDECOM(2,2)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               I/O, vectors       ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TDECOM(1,3)
         WTIM = TDECOM(2,3)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               vector subtraction ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         IF (PARALG) THEN
            CTIM = TDECOM(1,4)
            WTIM = TDECOM(2,4)
            IF (DOCPCT) CPCT = CTIM*CFAC
            IF (DOWPCT) WPCT = WTIM*WFAC
            WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &      '               qualified CD       ',
     &      CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         END IF
         CTIM = TMISC(1,1)
         WTIM = TMISC(2,1)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   'Misc.          qualification      ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TMISC(1,2)
         WTIM = TMISC(2,2)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               red. set write     ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         CTIM = TMISC(1,3)
         WTIM = TMISC(2,3)
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               info write         ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         IF (PARALG) THEN
            CTIM = TMISC(1,4)
            WTIM = TMISC(2,4)
            IF (DOCPCT) CPCT = CTIM*CFAC
            IF (DOWPCT) WPCT = WTIM*WFAC
            WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &      '               diagonal sync      ',
     &      CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
            CTIM = TMISC(1,5)
            WTIM = TMISC(2,5)
            IF (DOCPCT) CPCT = CTIM*CFAC
            IF (DOWPCT) WPCT = WTIM*WFAC
            WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &      '               vector count sync  ',
     &      CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         END IF
         CTIM = CMISC
         WTIM = WMISC
         IF (DOCPCT) CPCT = CTIM*CFAC
         IF (DOWPCT) WPCT = WTIM*WFAC
         WRITE(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)')
     &   '               etc.               ',
     &   CTIM/6.0D1,CPCT,WTIM/6.0D1,WPCT
         WRITE(LUPRI,'(A,A)')
     &   '---------------------------------------------------',
     &   '----------------------'
         WRITE(LUPRI,'(A,1X,F10.2,10X,F10.2)')
     &   'Total:                            ',
     &   TDECDRV(1)/6.0D1,TDECDRV(2)/6.0D1
         WRITE(LUPRI,'(A,A)')
     &   '---------------------------------------------------',
     &   '----------------------'
         WRITE(LUPRI,'(A,I12)')
     &   'Total #system calls for vector read  :',NSYS_CALL
         IF (.NOT. CHO_SSCREEN) THEN
            WRITE(LUPRI,'(A,I12)')
     &      'Total #DGEMM  calls for vector subtr.:',NDGM_CALL
         END IF

      END IF

C     Screening statistics from vector subtractions.
C     ----------------------------------------------

      IF (CHO_SSCREEN) THEN
         CALL CHO_HEAD('Screening Statistics from Vector Subtraction',
     &                 '-',80,LUPRI)
         WRITE(LUPRI,'(/,A,12X,A)')
     &   'Norm used for diagonals      : ',SSNORM
         WRITE(LUPRI,'(A,1P,D15.4)')
     &   'Screening threshold          : ',SSTAU
         WRITE(LUPRI,'(A,1P,D15.4)')
     &   'Maximum possible #DGEMV calls: ',SUBSCRSTAT(1)
         WRITE(LUPRI,'(A,1P,D15.4)')
     &   'Actual #DGEMV calls          : ',SUBSCRSTAT(2)
         IF (SUBSCRSTAT(1) .GT. 0.0D0) THEN
            SSCRPCT = 1.0D2*(SUBSCRSTAT(1)-SUBSCRSTAT(2))/SUBSCRSTAT(1)
            WRITE(LUPRI,'(A,8X,F7.2,A)')
     &      'Screening percent            : ',SSCRPCT,'%'
         END IF
      END IF

C     Statistics for shell quadruples spanned by each reduced set.
C     (This is for test purposes.)
C     ------------------------------------------------------------

      IF (CHO_TSTSCREEN .AND.
     &    .NOT.(CHO_FAKE_PAR.AND.NPROCS.GT.1.AND.Is_Real_Par())) THEN

         CHO_SSCREEN_SAVE = CHO_SSCREEN
         CHO_SSCREEN = .TRUE. ! to avoid error termination

         IF (NTAU .LT. 1) THEN
            WRITE(LUPRI,*) SECNAM,': screening test requested, but ',
     &                     'NTAU is non-positive!'
            WRITE(LUPRI,*) SECNAM,': test is skipped!'
            GO TO 199 ! skip
         END IF
         TAU(1) = THRCOM
         DO ITAU = 2,NTAU
            TAU(ITAU) = TAU(ITAU-1)*DTAU
         END DO

         CALL CHO_HEAD('RS Screening Statistics (TEST)','-',80,LUPRI)
         WRITE(LUPRI,'(/,A,A)')
     &   'Norm used for diagonal shell pairs: ',SSNORM
         CALL CHO_SUBSCR_INIT()

         ILOC = 3

         DO ISYM = 1,NSYM
            IF (NUMCHO(ISYM) .GT. 0) THEN

               LRED = INFVEC(NUMCHO(ISYM),2,ISYM)
               DO IRED = 1,LRED

                  IVEC1 = 0
                  NVEC  = 0
                  CALL CHO_X_NVECRS(IRED,ISYM,IVEC1,NVEC)
                  IF (NVEC.GT.0 .AND. NDIMRS(ISYM,IRED).GT.0) THEN

                     CALL CHO_MEM('TstS Max','GETM','REAL',KRDVT,LRDVT)
                     NUMVEC = MIN(LRDVT/NDIMRS(ISYM,IRED),NVEC)
                     IF (NUMVEC .LT. 1) THEN
                        CALL CHO_QUIT(
     &                  'Insufficient memory for TstScreen in '//SECNAM,
     &                  104)
                     END IF
                     NBATCH = (NVEC - 1)/NUMVEC + 1

                     DO IBATCH = 1,NBATCH

                        IF (IBATCH .EQ. NBATCH) THEN
                           NUMV = NVEC - NUMVEC*(NBATCH-1)
                        ELSE
                           NUMV = NUMVEC
                        END IF
                        LRDVEC = NDIMRS(ISYM,IRED)*NUMV
                        CALL CHO_MEM('TstScreen','ALLO','REAL',
     &                               KRDVEC,LRDVEC)

                        JVEC1 = IVEC1 + NUMVEC*(IBATCH-1)
                        JVEC2 = JVEC1 + NUMV - 1
                        JNUM  = 0
                        MUSD  = 0
                        CALL CHO_VECRD(WORK(KRDVEC),LRDVEC,JVEC1,JVEC2,
     &                                 ISYM,JNUM,IREDC,MUSD)
                        IF (JNUM .NE. NUMV) THEN
                           CALL CHO_QUIT('Logical error in '//SECNAM,
     &                                   103)
                        END IF

                        IF (IREDC .NE. IRED) THEN
                           CALL CHO_X_SETRED(IRC,ILOC,IRED)
                           IF (IRC .NE. 0) THEN
                              WRITE(LUPRI,*) SECNAM,
     &                        ': CHO_X_SETRED returned ',IRC
                              CALL CHO_QUIT('Error in '//SECNAM,104)
                           END IF
                           IREDC = IRED
                        END IF

                        CALL CHO_SUBSCR_DIA(WORK(KRDVEC),NUMV,ISYM,ILOC,
     &                                      SSNORM)
                        XT = 0.0D0
                        CALL CHO_DZERO(XC,NTAU)
                        DO ISHAB = 1,NNSHL
                           IF (NNBSTRSH(ISYM,ISHAB,ILOC) .GT. 0) THEN
                              DO ISHGD = ISHAB,NNSHL
                                 IF (NNBSTRSH(ISYM,ISHGD,ILOC) .GT. 0)
     &                           THEN
                                    TST = DSPNM(ISHAB)*DSPNM(ISHGD)
                                    TST = SQRT(TST)
                                    XT  = XT + 1.0D0
                                    DO ITAU = 1,NTAU
                                       IF (TST .GT. TAU(ITAU)) THEN
                                          XC(ITAU) = XC(ITAU) + 1.0D0
                                       END IF
                                    END DO
                                 END IF
                              END DO
                           END IF
                        END DO

                        WRITE(LUPRI,'(/,1X,A,I6,A,I2,A)')
     &                  '*** Statistics for reduced set',IRED,
     &                  '    Symmetry',ISYM,' ***'
                        WRITE(LUPRI,'(1X,A,I6,A,I6)')
     &                  '    Batch no.',IBATCH,' of',NBATCH
                        WRITE(LUPRI,'(1X,A,9X,I6)')
     &                  '       No. vectors     : ',NUMV
                        WRITE(LUPRI,'(1X,A,9X,I6,9X,I6)')
     &                  '       Vector range    : ',JVEC1,JVEC2
                        WRITE(LUPRI,'(1X,A,1P,D15.7)')
     &                  '       Shell quadruples: ',XT
                        IF (XT .LT. 1.0D0) THEN
                           CALL CHO_QUIT('XT non-positive in '//SECNAM,
     &                                   103)
                        END IF
                        FAC = 1.0D2/XT
                        DO ITAU = 1,NTAU
                           PCT = FAC*(XT-XC(ITAU))
                           WRITE(LUPRI,'(1X,A,1P,D15.7,A,D15.7,A)')
     &                     '       Threshold: ',TAU(ITAU),
     &                     '  Screening percent: ',PCT,'%'
                        END DO
                        CALL CHO_FLUSH(LUPRI)

                        CALL CHO_MEM('TstScreen','FREE','REAL',
     &                               KRDVEC,LRDVEC)

                     END DO

                  END IF

               END DO

            END IF
         END DO

         CALL CHO_SUBSCR_FINAL()

  199    CONTINUE ! skip point

         CHO_SSCREEN = CHO_SSCREEN_SAVE

      END IF


      END
