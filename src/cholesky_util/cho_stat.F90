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

subroutine CHO_STAT()
!
! Purpose: print statistics from decomposition.

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Para_Info, only: Is_Real_Par, nProcs
use Cholesky, only: CHO_DECALG, CHO_FAKE_PAR, CHO_INTCHK, CHO_REORD, Cho_SScreen, CHO_TSTSCREEN, DID_DECDRV, DSPNm, DSPNM, &
                    INF_INIT, INF_TIMING, InfVec, IntMap, IPRINT, LuPri, NBAS, NBAST, nDGM_call, nDimRS, nDimRS, nnBstR, nnBstRSh, &
                    nnBstRSh, nnBstRT, nnShl, nShell, nSym, nSys_call, NumCho, NumChT, RstCho, RstDia, SSNorm, SSTau, SubScrStat, &
                    TDECDRV, TDECOM, ThrCom, TIMSEC, TINTEG, TMISC, XnPass
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: NTAU = 5
integer(kind=iwp) :: IBATCH, ICAL, IHC, IHW, ILOC, IMC, IMW, IPRSAV, IRC, IRED, IREDC, ISHAB, ISHGD, ISHLAB, ISYM, ITAU, IVEC1, &
                     JNUM, JSYM, JVEC1, JVEC2, KSYM, LRDVEC, LRDVT, LRED, MAXCAL, MUSD, N, NBATCH, NCAL, NCALL, NN, NREP, NTOT, &
                     NUMV, NUMVEC, NVEC
real(kind=wp) :: CFAC, CMISC, CPCT, CTIM, FAC, PCT, SAV1, SAV2, SCC, SCW, SSCRPCT, TAU(NTAU), TCCHA, TCCHD, TCDEC, TCDIA, TCDIS, &
                 TCFIN, TCINI, TCREO, TST, TWCHA, TWCHD, TWDEC, TWDIA, TWDIS, TWFIN, TWINI, TWREO, VCSTOR(8), VCTOT, WFAC, WMISC, &
                 WPCT, WTIM, X, X1DIM, X1RED, X1TOT, X2TOT, XC(NTAU), XCAL, XGB, XREP, XT, XX, XXBST(8), XXSHL, XXTOT, XXX, Y, YY, &
                 YYY
logical(kind=iwp) :: CHO_SSCREEN_SAVE, DOCPCT, DOWPCT, PARALG
character(len=25) :: STRING
character(len=2) :: UNT
real(kind=wp), allocatable :: KRDVEC(:)
real(kind=wp), parameter :: DTAU = 1.0e-1_wp, MNT = 60.0_wp
character(len=*), parameter :: SECNAM = 'CHO_STAT'

PARALG = (CHO_DECALG == 4) .or. (CHO_DECALG == 5) .or. (CHO_DECALG == 6)

! Overall header.
! ---------------

call CHO_HEAD('Cholesky Decomposition Statistics','=',80,LUPRI)
if (RSTDIA) write(LUPRI,'(/,A)') 'Calculation restarted from diagonal on disk'
if (RSTCHO) then
  if (RSTDIA) then
    write(LUPRI,'(A)') 'Calculation restarted from Cholesky vectors on disk'
  else
    write(LUPRI,'(/,A)') 'Calculation restarted from Cholesky vectors on disk'
  end if
end if

! Configuration.
! --------------

call CHO_HEAD('Configuration','-',80,LUPRI)
IPRSAV = IPRINT
IPRINT = INF_INIT+1
call CHO_PRTHEAD(.true.)
IPRINT = IPRSAV

! Vector statistics.
! ------------------

call CHO_HEAD('Vector statistics','-',80,LUPRI)
write(LUPRI,'(/,A,/,A)') 'Sym.        N      Full      Mmax         M    M/Full    M/Mmax       M/N', &
                         '--------------------------------------------------------------------------'
NTOT = 0
do ISYM=1,NSYM
  NTOT = NTOT+NUMCHO(ISYM)
  NN = 0
  XXBST(ISYM) = Zero
  do JSYM=1,NSYM
    KSYM = MUL(JSYM,ISYM)
    if (JSYM > KSYM) then
      NN = NN+NBAS(JSYM)*NBAS(KSYM)
      XXBST(ISYM) = XXBST(ISYM)+real(NBAS(JSYM),kind=wp)*real(NBAS(KSYM),kind=wp)
    else if (JSYM == KSYM) then
      NN = NN+nTri_Elem(NBAS(JSYM))
      XXBST(ISYM) = XXBST(ISYM)+real(NBAS(JSYM),kind=wp)*(real(NBAS(JSYM),kind=wp)+One)*Half
    end if
  end do
  if (NN > 0) then
    YYY = real(NUMCHO(ISYM),kind=wp)/XXBST(ISYM)
  else
    YYY = 9.0e9_wp
  end if
  if (NNBSTR(ISYM,1) > 0) then
    XX = real(NNBSTR(ISYM,1),kind=wp)
    YY = real(NUMCHO(ISYM),kind=wp)/XX
  else
    YY = 9.0e9_wp
  end if
  if (NBAS(ISYM) /= 0) then
    X = real(NBAS(ISYM),kind=wp)
    Y = real(NUMCHO(ISYM),kind=wp)/X
  else
    Y = 9.0e9_wp
  end if
  write(LUPRI,'(I3,4(1X,I9),3(1X,F9.4))') ISYM,NBAS(ISYM),NN,NNBSTR(ISYM,1),NUMCHO(ISYM),YYY,YY,Y
end do
write(LUPRI,'(A)') '--------------------------------------------------------------------------'
NN = nTri_Elem(NBAST)
if (NN > 0) then
  XXX = real(NN,kind=wp)
  YYY = real(NUMCHT,kind=wp)/XXX
else
  YYY = 9.0e9_wp
end if
if (NNBSTRT(1) > 0) then
  XX = real(NNBSTRT(1),kind=wp)
  YY = real(NUMCHT,kind=wp)/XX
else
  YY = 9.0e9_wp
end if
if (NBAST /= 0) then
  X = real(NBAST,kind=wp)
  Y = real(NUMCHT,kind=wp)/X
else
  Y = 9.0e9_wp
end if
write(LUPRI,'(3X,4(1X,I9),3(1X,F9.4))') NBAST,NN,NNBSTRT(1),NUMCHT,YYY,YY,Y
write(LUPRI,'(A)') '--------------------------------------------------------------------------'
if (NTOT /= NUMCHT) write(LUPRI,'(A)') 'WARNING: total number of vectors is wrong!!!'

call CHO_GETSTOR(VCSTOR)

write(LUPRI,'(/,A,/,A,/,A)') '                       %Saving relative to','Sym.     Storage      1st Red. Set     Full', &
                             '---------------------------------------------'
VCTOT = Zero
X1TOT = Zero
X2TOT = Zero
XXTOT = Zero
do ISYM=1,NSYM
  X1DIM = real(NNBSTR(ISYM,1),kind=wp)
  X1RED = X1DIM*real(NUMCHO(ISYM),kind=wp)
  if (X1RED > Zero) then
    SAV1 = 1.0e2_wp*(X1RED-VCSTOR(ISYM))/X1RED
  else
    SAV1 = 9.0e9_wp
  end if
  XX = XXBST(ISYM)*real(NUMCHO(ISYM),kind=wp)
  if (XX > Zero) then
    SAV2 = 1.0e2_wp*(XX-VCSTOR(ISYM))/XX
  else
    SAV2 = 9.0e9_wp
  end if
  XGB = VCSTOR(ISYM)*Eight/1.024e3_wp
  UNT = 'kb'
  if (XGB > 1.0e3_wp) then
    XGB = XGB/1.024e3_wp
    UNT = 'Mb'
    if (XGB > 1.0e3_wp) then
      XGB = XGB/1.024e3_wp
      UNT = 'Gb'
      if (XGB > 1.0e3) then
        XGB = XGB/1.024e3_wp
        UNT = 'Tb'
      end if
    end if
  end if
  write(LUPRI,'(I2,4X,F10.3,1X,A,4X,F9.4,4X,F9.4)') ISYM,XGB,UNT,SAV1,SAV2
  VCTOT = VCTOT+VCSTOR(ISYM)
  X1TOT = X1TOT+X1RED
  X2TOT = X2TOT+X1DIM*(X1DIM+One)*Half
  XXTOT = XXTOT+XX
end do
write(LUPRI,'(A)') '---------------------------------------------'
if (X1TOT > Zero) then
  SAV1 = 1.0e2_wp*(X1TOT-VCTOT)/X1TOT
else
  SAV1 = 9.0e9_wp
end if
if (XXTOT > Zero) then
  SAV2 = 1.0e2_wp*(XXTOT-VCTOT)/XXTOT
else
  SAV2 = 9.0e9_wp
end if
XGB = VCTOT*Eight/1.024e3_wp
UNT = 'kb'
if (XGB > 1.0e3_wp) then
  XGB = XGB/1.024e3_wp
  UNT = 'Mb'
  if (XGB > 1.0e3_wp) then
    XGB = XGB/1.024e3_wp
    UNT = 'Gb'
    if (XGB > 1.0e3_wp) then
      XGB = XGB/1.024e3_wp
      UNT = 'Tb'
    end if
  end if
end if
write(LUPRI,'(A6,F10.3,1X,A,4X,F9.4,4X,F9.4)') 'Total:',XGB,UNT,SAV1,SAV2
write(LUPRI,'(A)') '---------------------------------------------'
PCT = 1.0e2_wp*VCTOT/X2TOT
write(LUPRI,'(A,F11.6,A)') 'Total storage corresponds to',PCT,'% of the 1st reduced set integral matrix.'

call CHO_STAT_PARENTDIAG()

! Integral statistics.
! --------------------

MAXCAL = 0
NCAL = 0
NREP = 0
NCALL = 0
do ISHLAB=1,NNSHL
  if (INTMAP(ISHLAB) > 0) then
    MAXCAL = max(MAXCAL,INTMAP(ISHLAB))
    NCAL = NCAL+1
    NCALL = NCALL+INTMAP(ISHLAB)
    if (INTMAP(ISHLAB) > 1) NREP = NREP+1
  end if
end do

XXSHL = real(NNSHL,kind=wp)
XCAL = real(NCAL,kind=wp)
XREP = real(NREP,kind=wp)

call CHO_HEAD('Integral statistics','-',80,LUPRI)
write(LUPRI,'(/,A,I10)') '#Shells                  :',NSHELL
write(LUPRI,'(A,I10)') '#Shell Pair Distributions:',NNSHL
write(LUPRI,'(A,I10)') '#Integral passes         :',XNPASS
write(LUPRI,'(A,I10)') '#Calls to integral prog. :',NCALL
write(LUPRI,'(A,I10,A,F8.3,A)') '#Shell Pairs Calculated  :',NCAL,' (',XCAL*1.0e2_wp/XXSHL,' % of total)'
write(LUPRI,'(A,I10,A,F8.3,A)') '#Shell Pairs Repeated    :',NREP,' (',XREP*1.0e2_wp/XCAL,' % of calculated)'

write(LUPRI,'(/,A,/,A)') '#Calculations     #Shell Pairs   Percentage','-------------------------------------------'
do ICAL=1,MAXCAL
  N = 0
  do ISHLAB=1,NNSHL
    if (INTMAP(ISHLAB) == ICAL) N = N+1
  end do
  if (N > 0) then
    X = real(N,kind=wp)
    write(LUPRI,'(I12,6X,I12,5X,F8.3)') ICAL,N,X*1.0e2_wp/XXSHL
  end if
end do
write(LUPRI,'(A)') '-------------------------------------------'

! Section timings.
! ----------------

if (IPRINT >= INF_TIMING) then

  TCINI = TIMSEC(2,1)-TIMSEC(1,1)
  TWINI = TIMSEC(4,1)-TIMSEC(3,1)
  TCDIA = TIMSEC(2,2)-TIMSEC(1,2)
  TWDIA = TIMSEC(4,2)-TIMSEC(3,2)
  TCDEC = TIMSEC(2,3)-TIMSEC(1,3)
  TWDEC = TIMSEC(4,3)-TIMSEC(3,3)
  TCCHD = TIMSEC(2,4)-TIMSEC(1,4)
  TWCHD = TIMSEC(4,4)-TIMSEC(3,4)
  TCCHA = TIMSEC(2,5)-TIMSEC(1,5)
  TWCHA = TIMSEC(4,5)-TIMSEC(3,5)
  TCREO = TIMSEC(2,6)-TIMSEC(1,6)
  TWREO = TIMSEC(4,6)-TIMSEC(3,6)
  TCDIS = TIMSEC(2,7)-TIMSEC(1,7)
  TWDIS = TIMSEC(4,7)-TIMSEC(3,7)
  TCFIN = TIMSEC(2,8)-TIMSEC(1,8)
  TWFIN = TIMSEC(4,8)-TIMSEC(3,8)

  call CHO_HEAD('Section timings','-',80,LUPRI)
  write(LUPRI,'(/,A,/,A,/,A)') '                                 CPU time          Wall time', &
                               'Section                      hours min. sec.    hours min. sec.', &
                               '---------------------------------------------------------------'
  STRING = 'Initialization           '
  call CHO_CNVTIM(TCINI,IHC,IMC,SCC)
  call CHO_CNVTIM(TWINI,IHW,IMW,SCW)
  write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  STRING = 'Diagonal setup           '
  call CHO_CNVTIM(TCDIA,IHC,IMC,SCC)
  call CHO_CNVTIM(TWDIA,IHW,IMW,SCW)
  write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  STRING = 'Cholesky decomposition   '
  call CHO_CNVTIM(TCDEC,IHC,IMC,SCC)
  call CHO_CNVTIM(TWDEC,IHW,IMW,SCW)
  write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  STRING = 'Diagonal check           '
  call CHO_CNVTIM(TCCHD,IHC,IMC,SCC)
  call CHO_CNVTIM(TWCHD,IHW,IMW,SCW)
  write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  if (CHO_INTCHK) then
    STRING = 'Integral check (debug)   '
    call CHO_CNVTIM(TCCHA,IHC,IMC,SCC)
    call CHO_CNVTIM(TWCHA,IHW,IMW,SCW)
    write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  end if
  if (CHO_REORD) then
    STRING = 'Vector reordering        '
    call CHO_CNVTIM(TCREO,IHC,IMC,SCC)
    call CHO_CNVTIM(TWREO,IHW,IMW,SCW)
    write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  end if
  if (CHO_FAKE_PAR .and. (NPROCS > 1) .and. Is_Real_Par()) then
    STRING = 'Vector distribution      '
    call CHO_CNVTIM(TCDIS,IHC,IMC,SCC)
    call CHO_CNVTIM(TWDIS,IHW,IMW,SCW)
    write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  end if
  STRING = 'Finalization             '
  call CHO_CNVTIM(TCFIN,IHC,IMC,SCC)
  call CHO_CNVTIM(TWFIN,IHW,IMW,SCW)
  write(LUPRI,'(A,1X,I8,2X,I2,1X,F5.1,1X,I8,2X,I2,1X,F5.1)') STRING,IHC,IMC,SCC,IHW,IMW,SCW
  write(LUPRI,'(A)') '---------------------------------------------------------------'

end if ! IPRINT > INF_TIMING

! Timing of decomposition driver.
! -------------------------------

if (DID_DECDRV) then

  CMISC = TDECDRV(1)-sum(TINTEG(1,:))-sum(TDECOM(1,:))-sum(TMISC(1,:))
  WMISC = TDECDRV(2)-sum(TINTEG(2,:))-sum(TDECOM(2,:))-sum(TMISC(2,:))
  CPCT = -9.0e9_wp
  CFAC = -9.0e9_wp
  WPCT = -9.0e9_wp
  WFAC = -9.0e9_wp
  if (TDECDRV(1) > Zero) then
    DOCPCT = .true.
    CFAC = 1.0e2_wp/TDECDRV(1)
  else
    DOCPCT = .false.
  end if
  if (TDECDRV(2) > Zero) then
    DOWPCT = .true.
    WFAC = 1.0e2_wp/TDECDRV(2)
  else
    DOWPCT = .false.
  end if

  call CHO_HEAD('Timing of decomposition driver','-',80,LUPRI)
  write(LUPRI,'(/,A,/,A)') 'Task           Component           CPU (min.)     %   Wall (min.)     %', &
                           '-------------------------------------------------------------------------'
  CTIM = TINTEG(1,1)
  WTIM = TINTEG(2,1)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') 'Integrals      calculation        ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TINTEG(1,2)
  WTIM = TINTEG(2,2)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               I/O, qualifieds    ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TDECOM(1,1)
  WTIM = TDECOM(2,1)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') 'Decomposition  I/O, qualifieds    ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TDECOM(1,2)
  WTIM = TDECOM(2,2)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               I/O, vectors       ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TDECOM(1,3)
  WTIM = TDECOM(2,3)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               vector subtraction ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  if (PARALG) then
    CTIM = TDECOM(1,4)
    WTIM = TDECOM(2,4)
    if (DOCPCT) CPCT = CTIM*CFAC
    if (DOWPCT) WPCT = WTIM*WFAC
    write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               qualified CD       ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  end if
  CTIM = TMISC(1,1)
  WTIM = TMISC(2,1)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') 'Misc.          qualification      ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TMISC(1,2)
  WTIM = TMISC(2,2)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               red. set write     ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  CTIM = TMISC(1,3)
  WTIM = TMISC(2,3)
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               info write         ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  if (PARALG) then
    CTIM = TMISC(1,4)
    WTIM = TMISC(2,4)
    if (DOCPCT) CPCT = CTIM*CFAC
    if (DOWPCT) WPCT = WTIM*WFAC
    write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               diagonal sync      ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
    CTIM = TMISC(1,5)
    WTIM = TMISC(2,5)
    if (DOCPCT) CPCT = CTIM*CFAC
    if (DOWPCT) WPCT = WTIM*WFAC
    write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               vector count sync  ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  end if
  CTIM = CMISC
  WTIM = WMISC
  if (DOCPCT) CPCT = CTIM*CFAC
  if (DOWPCT) WPCT = WTIM*WFAC
  write(LUPRI,'(A,1X,F10.2,1X,F7.2,2X,F10.2,1X,F7.2)') '               etc.               ',CTIM/MNT,CPCT,WTIM/MNT,WPCT
  write(LUPRI,'(A)') '-------------------------------------------------------------------------'
  write(LUPRI,'(A,1X,F10.2,10X,F10.2)') 'Total:                            ',TDECDRV(1)/MNT,TDECDRV(2)/MNT
  write(LUPRI,'(A)') '-------------------------------------------------------------------------'
  write(LUPRI,'(A,I12)') 'Total #system calls for vector read  :',NSYS_CALL
  if (.not. CHO_SSCREEN) write(LUPRI,'(A,I12)') 'Total #DGEMM  calls for vector subtr.:',NDGM_CALL

end if

! Screening statistics from vector subtractions.
! ----------------------------------------------

if (CHO_SSCREEN) then
  call CHO_HEAD('Screening Statistics from Vector Subtraction','-',80,LUPRI)
  write(LUPRI,'(/,A,12X,A)') 'Norm used for diagonals      : ',SSNORM
  write(LUPRI,'(A,1P,D15.4)') 'Screening threshold          : ',SSTAU
  write(LUPRI,'(A,1P,D15.4)') 'Maximum possible #DGEMV calls: ',SUBSCRSTAT(1)
  write(LUPRI,'(A,1P,D15.4)') 'Actual #DGEMV calls          : ',SUBSCRSTAT(2)
  if (SUBSCRSTAT(1) > Zero) then
    SSCRPCT = 1.0e2_wp*(SUBSCRSTAT(1)-SUBSCRSTAT(2))/SUBSCRSTAT(1)
    write(LUPRI,'(A,8X,F7.2,A)') 'Screening percent            : ',SSCRPCT,'%'
  end if
end if

! Statistics for shell quadruples spanned by each reduced set.
! (This is for test purposes.)
! ------------------------------------------------------------

if (CHO_TSTSCREEN .and. (.not. (CHO_FAKE_PAR .and. (NPROCS > 1) .and. Is_Real_Par()))) then

  CHO_SSCREEN_SAVE = CHO_SSCREEN
  CHO_SSCREEN = .true. ! to avoid error termination

  if (NTAU < 1) then
    write(LUPRI,*) SECNAM,': screening test requested, but NTAU is non-positive!'
    write(LUPRI,*) SECNAM,': test is skipped!'
  else
    TAU(1) = THRCOM
    do ITAU=2,NTAU
      TAU(ITAU) = TAU(ITAU-1)*DTAU
    end do

    call CHO_HEAD('RS Screening Statistics (TEST)','-',80,LUPRI)
    write(LUPRI,'(/,A,A)') 'Norm used for diagonal shell pairs: ',SSNORM
    call CHO_SUBSCR_INIT()

    ILOC = 3

    do ISYM=1,NSYM
      if (NUMCHO(ISYM) > 0) then

        LRED = INFVEC(NUMCHO(ISYM),2,ISYM)
        do IRED=1,LRED

          IVEC1 = 0
          NVEC = 0
          call CHO_X_NVECRS(IRED,ISYM,IVEC1,NVEC)
          if ((NVEC > 0) .and. (NDIMRS(ISYM,IRED) > 0)) then

            call mma_maxDBLE(LRDVT)
            NUMVEC = min(LRDVT/NDIMRS(ISYM,IRED),NVEC)
            if (NUMVEC < 1) call CHO_QUIT('Insufficient memory for TstScreen in '//SECNAM,104)
            NBATCH = (NVEC-1)/NUMVEC+1

            do IBATCH=1,NBATCH

              if (IBATCH == NBATCH) then
                NUMV = NVEC-NUMVEC*(NBATCH-1)
              else
                NUMV = NUMVEC
              end if
              LRDVEC = NDIMRS(ISYM,IRED)*NUMV
              call mma_allocate(KRDVEC,LRDVEC,Label='KRDVEC')

              JVEC1 = IVEC1+NUMVEC*(IBATCH-1)
              JVEC2 = JVEC1+NUMV-1
              JNUM = 0
              MUSD = 0
              call CHO_VECRD(KRDVEC,LRDVEC,JVEC1,JVEC2,ISYM,JNUM,IREDC,MUSD)
              if (JNUM /= NUMV) call CHO_QUIT('Logical error in '//SECNAM,103)

              if (IREDC /= IRED) then
                call CHO_X_SETRED(IRC,ILOC,IRED)
                if (IRC /= 0) then
                  write(LUPRI,*) SECNAM,': CHO_X_SETRED returned ',IRC
                  call CHO_QUIT('Error in '//SECNAM,104)
                end if
                IREDC = IRED
              end if

              call CHO_SUBSCR_DIA(KRDVEC,NUMV,ISYM,ILOC,SSNORM)
              XT = Zero
              XC(:) = Zero
              do ISHAB=1,NNSHL
                if (NNBSTRSH(ISYM,ISHAB,ILOC) > 0) then
                  do ISHGD=ISHAB,NNSHL
                    if (NNBSTRSH(ISYM,ISHGD,ILOC) > 0) then
                      TST = DSPNM(ISHAB)*DSPNM(ISHGD)
                      TST = sqrt(TST)
                      XT = XT+One
                      do ITAU=1,NTAU
                        if (TST > TAU(ITAU)) XC(ITAU) = XC(ITAU)+One
                      end do
                    end if
                  end do
                end if
              end do

              write(LUPRI,'(/,1X,A,I6,A,I2,A)') '*** Statistics for reduced set',IRED,'    Symmetry',ISYM,' ***'
              write(LUPRI,'(1X,A,I6,A,I6)') '    Batch no.',IBATCH,' of',NBATCH
              write(LUPRI,'(1X,A,9X,I6)') '       No. vectors     : ',NUMV
              write(LUPRI,'(1X,A,9X,I6,9X,I6)') '       Vector range    : ',JVEC1,JVEC2
              write(LUPRI,'(1X,A,1P,D15.7)') '       Shell quadruples: ',XT
              if (XT < One) call CHO_QUIT('XT non-positive in '//SECNAM,103)
              FAC = 1.0e2_wp/XT
              do ITAU=1,NTAU
                PCT = FAC*(XT-XC(ITAU))
                write(LUPRI,'(1X,A,1P,D15.7,A,D15.7,A)') '       Threshold: ',TAU(ITAU),'  Screening percent: ',PCT,'%'
              end do
              call XFLUSH(LUPRI)

              call mma_deallocate(KRDVEC)

            end do

          end if

        end do

      end if
    end do

    call CHO_SUBSCR_FINAL()
  end if ! skip point

  CHO_SSCREEN = CHO_SSCREEN_SAVE

end if

end subroutine CHO_STAT
