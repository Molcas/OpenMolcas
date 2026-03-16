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

subroutine READIN_RASSI()

use rassi_aux, only: ipglob
use rassi_global_arrays, only: HAM, ESHFT, HDIAG, JBNUM, LROOT
use frenkel_global_vars, only: excl, iTyp, valst, corest, nesta, nestb, nestla, nestlb, doexch, DoExcitonics, DoCoul, labA, labB, &
                               rixs
use kVectors, only: e_Vector, k_Vector, nk_Vector
use Lebedev_quadrature, only: available_table, rule_max
#ifdef _DMRG_
use rasscf_global, only: doDMRG
use qcmaquis_interface_cfg
#endif
use Fock_util_global, only: Deco, Estimate, PseudoChoMOs, Update
use Cholesky, only: timings
use stdalloc, only: mma_allocate
use cntrl, only: SONTO, SONTOSTATES, SONAT, SONATNSTATE, SODIAG, SODIAGNSTATE
use spool, only: Spoolinp, Close_LuSpool
use Cntrl, only: QDPT2SC, QDPT2EV, SECOND_TIME, DOGSOR, PRSXY, PRORB, PRTRA, PRCI, BINA, NATO, NBINA, NRNATO, RFPERT, IFTRD1, &
                 NSOPR, NPROP, PRDIPVEC, TDIPMIN, NJOB, CIH5, CIThr, IFHAM, IFSO, IFNTO, SOThr_Prt, nSOThr_Prt, nState, IfHEXT, &
                 IfHEff, IfHCOM, IFEJOB, IfHDia, IfShft, ToFile, IfJ2, IfJZ, IFGCAL, EPraThr, IFACALSD, IFACALFC, IFACALSDON, &
                 IFACALPSO, IFATCALSA, IFGTSHSA, MULTIP, IFVANVLECK, TMINS, TMAXS, NTS, IFSONCINI, TMINP, TMAXP, NTP, IFSONCIFC, &
                 TMINF, TMAXF, NTF, NOSO, IFCURD, IFARGU, IFXCAL, NBSTEP, BSTART, BINCRE, BANGRES, NTSTEP, TSTART, TINCRE, IFMCAL, &
                 PRXVR, PRXVE, PRXVS, PRMER, PRMEE, PRMES, HOP, TRACK, NOHAM, ONLY_OVERLAPS, IFDCPL, IFTRD2, IFTDM, DQVD, ALPHZ, &
                 BETAE, DIPR, OSTHR_DIPR, QIPR, OSTHR_QIPR, QIALL, RSPR, RSThr, DOCD, DYSO, DYSEXPORT, DYSEXPSO, TDYS, OCAN, DCHS, &
                 DCHO, DO_TMOM, TMGR_Thrs, PRRAW, PRWEIGHT, TOLERANCE, REDUCELOOP, LOOPDIVIDE, LOOPMAX, l_Eff, Do_SK, Do_Pol, &
                 RHODYN, MXJOB, JBNAME, SOPRNM, PNAME, PRDIPCOM, EPrThr, LPRPR, lHami, IfACAL, IFACALFCON, IFACALFCSDON, &
                 IFGTCALSA, DYSEXPSF, ISTAT, MXPROP, NSTAT, IBINA, ISOCMP, ICOMP, OCAA
use cntrl, only: ALGO, Nscreen, dmpk
use rassi_data, only: CHFRACMEM
use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit none
character(len=80) LINE
integer, parameter :: MXPLST = 50
character(len=8) TRYNAME
real*8 tmp
logical lExists
integer I, J, ISTATE, JSTATE, IJOB, LINENR
integer LuIn
integer NFLS
character(len=7) :: input_id

call SpoolInp(LuIn)

! Default settings for Cholesky
Algo = 2
Nscreen = 10
dmpk = 0.1_wp
timings = .false.
Estimate = .false.
Update = .true.
Deco = .true.
PseudoChoMOs = .false.
#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Zero
#endif

ITYP = 0
VALST = 1
COREST = 0
!> set some defaults for MPSSI
QDPT2SC = .true.
QDPT2EV = .false.
input_id = '&RASSI'
#ifdef _DMRG_
!> make sure that we read checkpoint names from xxx.h5 files, for example: rasscf.h5, nevpt2.h5, caspt2.h5, ...
if (doDMRG) input_id = '&MPSSI'
#endif

! Defaults for SI-PDFT runs:
Second_time = .false.
DoGSOR = .false.

! Find beginning of input:
50 read(LuIn,'(A72)',end=998) LINE
call NORMAL(LINE)
if (LINE(1:7) /= input_id) goto 50
LINENR = 0
100 read(LuIn,'(A72)',end=998) LINE
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') goto 100
if (LINE == '') goto 100
if (LINE(1:4) == 'END ') goto 200
! ------------------------------------------
if (LINE(1:4) == 'TEST') then
  PRSXY = .true.
  PRORB = .true.
  PRTRA = .true.
  PRCI = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'SECO') then
  SECOND_TIME = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'GSOR') then
  DoGSOR = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'BINA') then
  BINA = .true.
  NATO = .true.
  read(LuIn,*,err=997) NBINA
  LINENR = LINENR+1
  read(LuIn,*,err=997) (IBINA(1,I),IBINA(2,I),I=1,NBINA)
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'EXTR') then
  if (IPGLOB > 0) then
    call WarningMessage(1,'Obsolete EXTRACT keyword used.')
    write(u6,*) ' Please remove it from the input.'
  end if
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'NATO') then
  NATO = .true.
  read(LuIn,*,err=997) NRNATO
  LINENR = LINENR+1
  goto 100
end if
!-------------------------------------------
if (LINE(1:4) == 'RFPE') then
  RFPERT = .true.
  goto 100
end if
! -- FA 2005 start--------------------------
!-------------------------------------------
! --- Cholesky with default settings
if (LINE(1:4) == 'CHOL') then
  call Cho_rassi_rdInp(.true.,LuIn)
  goto 100
end if
! --- Cholesky with customized settings
if (LINE(1:4) == 'CHOI') then
  call Cho_rassi_rdInp(.false.,LuIn)
  goto 100
end if
! -- FA 2005 end----------------------------
if (LINE(1:4) == 'EXAL') then
  EXCL = .true.
  NESTA = 0
  read(LuIn,*,err=997) NESTA
  call mma_allocate(NESTLA,NESTA)
  LINENR = LINENR+1
  read(LuIn,*,err=997) (NESTLA(I),I=1,NESTA)
  !write(u6,*) 'list of excit. initial states'
  !write(u6,*) 'mon A nr states', NESTA
  !write(u6,'(I2)') (NESTLA(I),I=1,NESTA)
  goto 100
end if

if (LINE(1:4) == 'EXBL') then
  EXCL = .true.
  NESTB = 0
  read(LuIn,*,err=997) NESTB
  call mma_allocate(NESTLB,NESTB)
  LINENR = LINENR+1
  read(LuIn,*,err=997) (NESTLB(I),I=1,NESTB)
  !write(u6,*) 'list of excit. initial states'
  !write(u6,*) 'mon B nr states', NESTB
  !write(u6,'(I2)') (NESTLB(I),I=1,NESTB)
  goto 100
end if

if (LINE(1:4) == 'EXCI') then
  DoExcitonics = .true.
  goto 100
end if
if (LINE(1:4) == 'KCOU') then
  DoExch = .true.
  goto 100
end if
if (LINE(1:4) == 'MONA') then
  DoCoul = .true.
  labA = .true.
  if (iTyp == 2) write(u6,*) ' Warning: switching to monomer-A.'
  iTyp = 1
  if (.not. IFTRD1) then
    IFTRD1 = .true.
    LINENR = LINENR+1
    write(u6,*) ' TRD1 activated by MONA.'
  end if
  goto 100
end if
if (LINE(1:4) == 'MONB') then
  labB = .true.
  DoCoul = .true.
  if (iTyp == 1) write(u6,*) ' Warning: switching to monomer-B.'
  iTyp = 2
  if (.not. IFTRD1) then
    IFTRD1 = .true.
    LINENR = LINENR+1
    write(u6,*) ' TRD1 activated by MONB.'
  end if
  goto 100
end if
if (LINE(1:4) == 'RIXS') then
  RIXS = .true.
  read(LuIn,*,err=997) valst,corest
  LINENR = LINENR+1
  goto 100
end if
! --- FA 2016 end---------------------------
if (LINE(1:4) == 'SOPR') then
  read(LuIn,*,err=997) NSOPR,(SOPRNM(I),ISOCMP(I),I=1,min(MXPROP,NSOPR))
  LINENR = LINENR+1
  do I=1,min(MXPROP,NSOPR)
    call UPCASE(SOPRNM(I))
    if (SOPRNM(I)(1:5) == 'MLTPL') then
      if (SOPRNM(I)(7:8) == '  ') then
        SOPRNM(I) = 'MLTPL  '//SOPRNM(I)(6:6)
      else if (SOPRNM(I)(8:8) == ' ') then
        SOPRNM(I) = 'MLTPL '//SOPRNM(I)(6:7)
      end if
    end if
  end do
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'PROP') then
  read(LuIn,*,err=997) NPROP,(PNAME(I),ICOMP(I),I=1,min(MXPROP,NPROP))
  LINENR = LINENR+1
  do I=1,min(MXPROP,NPROP)
    call UPCASE(PNAME(I))
    if (PNAME(I)(1:5) == 'MLTPL') then
      if (PNAME(I)(7:8) == '  ') then
        PNAME(I) = 'MLTPL  '//PNAME(I)(6:6)
      else if (PNAME(I)(8:8) == ' ') then
        PNAME(I) = 'MLTPL '//PNAME(I)(6:7)
      end if
    end if
  end do
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'OVER') then
  PRSXY = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'TRDI') then
  ! Print transition dipole vectors
  PRDIPVEC = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'TDMN') then
  ! Print transition dipole vectors
  PRDIPVEC = .true.
  read(LuIn,*,err=997) TDIPMIN
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'TRDC') then
  ! Print COMPLEX transition dipole vectors
  PRDIPCOM = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'ORBI') then
  PRORB = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'CIPR') then
  PRCI = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'CIH5') then
  if (NJOB <= 2) then
    CIH5 = .true.
  else
    call WarningMessage(2,'CIH5 allows no more than 2 JOBIPHs')
    call abend()
  end if
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'THRS') then
  read(LuIn,*,err=997) CITHR
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'ONEL') then
  IFHAM = .false.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'ONEE') then
  IFHAM = .false.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'SPIN') then
  IFSO = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'NTOC') then
  IFNTO = .true.
  goto 100
end if
! ------------------------------------------
! PAM07 Added: Keyword for printing spin-orbit coupling matrix elements
! A threshold in reciprocal cm is entered.
if (LINE(1:4) == 'SOCO') then
  read(LuIn,*,err=997) SOTHR_PRT
  LINENR = LINENR+1
  if (SOTHR_PRT < Zero) SOTHR_PRT = Zero
  NSOTHR_PRT = 10000
  goto 100
end if
! ------------------------------------------
if ((LINE(1:4) == 'NROF') .or. (LINE(1:4) == 'NR O')) then
  NSTATE = 0
  read(LuIn,*,err=997) NJOB,TRYNAME
  call UpCase(TRYNAME)
  if (TRYNAME == 'ALL') then
    LINENR = LINENR+1
  else
    backspace(LuIn)
    read(LuIn,*,err=997) NJOB,(NSTAT(I),I=1,NJOB)
    do IJOB=1,NJOB
      NSTATE = NSTATE+NSTAT(IJOB)
    end do
    call mma_allocate(JBNUM,nState,Label='JBNUM')
    call mma_allocate(LROOT,nState,Label='LROOT')
    LINENR = LINENR+1
    NSTATE = 0
    do IJOB=1,NJOB
      ISTAT(IJOB) = NSTATE+1
      read(LuIn,*,err=997) (LROOT(NSTATE+J),J=1,NSTAT(IJOB))
      LINENR = LINENR+1
      do ISTATE=NSTATE+1,NSTATE+NSTAT(IJOB)
        JBNUM(ISTATE) = IJOB
      end do
      NSTATE = NSTATE+NSTAT(IJOB)
    end do
  end if
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'IPHN') then
  !read(LuIn,'(9(A7,1X))',err=997) (JBNAME(I),I=1,NJOB)
  !LINENR = LINENR+1
  do I=1,NJOB
    read(LuIn,*,err=997) JBNAME(I)
    LINENR = LINENR+1
    call UpCase(JBNAME(I))
  end do
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'FILE') then
  read(LuIn,*,err=997) NJOB
  LINENR = LINENR+1
  do I=1,NJOB
    read(LuIn,*,err=997) LINE
    LINENR = LINENR+1
    call FILEORB(LINE,JBNAME(I))
  end do
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'HEXT') then
  if (nState == 0) then
    call WarningMessage(2,'HEXT needs an explicit NROFJOB')
    call abend()
  end if
  IFHEXT = .true.
  call mma_allocate(HAM,nState,nState,Label='HAM')
  read(LuIn,*,err=997) ((HAM(ISTATE,JSTATE),JSTATE=1,ISTATE),ISTATE=1,NSTATE)
  do ISTATE=1,NSTATE-1
    do JSTATE=ISTATE+1,NSTATE
      HAM(ISTATE,JSTATE) = HAM(JSTATE,ISTATE)
    end do
  end do
  LINENR = LINENR+NSTATE
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'HEFF') then
  IFHEFF = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'HCOM') then
  IFHCOM = .true.
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'EJOB') then
  IFEJOB = .true.
  !   Leon: Is it really needed?
  !LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'HDIA') then
  if (nState == 0) then
    call WarningMessage(2,'HDIA needs an explicit NROFJOB')
    call abend()
  end if
  IFHDIA = .true.
  call mma_allocate(HDIAG,nState,Label='nState')
  read(LuIn,*,err=997) (HDIAG(ISTATE),ISTATE=1,NSTATE)
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'SHIF') then
  if (nState == 0) then
    call WarningMessage(2,'SHIFT needs an explicit NROFJOB')
    call abend()
  end if
  IFSHFT = .true.
  call mma_allocate(ESHFT,nState,Label='ESHFT')
  read(LuIn,*,err=997) (ESHFT(ISTATE),ISTATE=1,NSTATE)
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'TOFI') then
  ToFile = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'J-VA') then
  IFJ2 = 1
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'OMEG') then
  IFJZ = 1
  Linenr = Linenr+1
  goto 100
end if
!-SVC 2007-----------------------------------
if (Line(1:4) == 'EPRG') then
  IFGCAL = .true.
  read(LuIn,*,err=997) EPRTHR
  if (EPRTHR < Zero) EPRTHR = Zero
  Linenr = Linenr+1
  goto 100
end if

! tjd- BMII: Print out spin-orbit properties to files
if (Line(1:4) == 'PRPR') then
  write(u6,*) 'SPIN-ORBIT PROPERTY PRINT ON'
  LPRPR = .true.
  Linenr = Linenr+1
  goto 100
end if

! tjd- Yoni: Force an identity SO hamiltonian
if (LINE(1:4) == 'HAMI') then
  write(u6,*) 'Identity Hamiltonian turned on'
  LHAMI = .true.
  Linenr = Linenr+1
  goto 100
end if

! BP - Hyperfine calculations
if (Line(1:4) == 'EPRA') then
  !write(u6,*)"EPRA read"
  IFACAL = .true.
  read(LuIn,*,err=997) EPRATHR
  if (EPRATHR < Zero) EPRATHR = Zero
  Linenr = Linenr+1
  goto 100
end if
if (Line(1:4) == 'AFCO') then
  !write(u6,*) 'AFCO read'
  IFACALSD = .false.
  Linenr = Linenr+1
  goto 100
end if
if (Line(1:4) == 'ASDO') then
  !write(u6,*) 'ASDO read'
  IFACALFC = .false.
  Linenr = Linenr+1
  goto 100
end if
! Kamal Sharkas beg - PSO Hyperfine calculations
if (Line(1:4) == 'AFCC') then
  IFACALFCON = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'ASDC') then
  IFACALSDON = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'FCSD') then
  IFACALFCSDON = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'APSO') then
  IFACALPSO = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'GTSA') then
  IFGTCALSA = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'ATSA') then
  IFATCALSA = .true.
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'SHMP') then
  IFGTSHSA = .true.
  read(LuIn,*,err=997) MULTIP
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'VSUS') then
  IFVANVLECK = .true.
  read(LuIn,*,err=997) TMINS,TMAXS,NTS
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'NMRT') then
  IFSONCINI = .true.
  read(LuIn,*,err=997) TMINP,TMAXP,NTP
  Linenr = Linenr+1
  goto 100
end if

if (Line(1:4) == 'NMRF') then
  IFSONCIFC = .true.
  read(LuIn,*,err=997) TMINF,TMAXF,NTF
  Linenr = Linenr+1
  goto 100
end if

! Kamal Sharkas end - PSO Hyperfine calculations

! BP Natural orbitals options
if (Line(1:4) == 'SONO') then
  read(LuIn,*,err=997) SONATNSTATE
  call mma_allocate(SONAT,SONATNSTATE,Label='SONAT')
  Linenr = Linenr+1
  read(LuIn,*,err=997) (SONAT(I),I=1,SONATNSTATE)
  Linenr = Linenr+1
  goto 100
end if
if (Line(1:4) == 'SODI') then
  read(LuIn,*,err=997) SODIAGNSTATE
  call mma_allocate(SODIAG,SODIAGNSTATE,Label='SODIAG')
  Linenr = Linenr+1
  read(LuIn,*,err=997) (SODIAG(I),I=1,SODIAGNSTATE)
  Linenr = Linenr+1
  goto 100
end if
if (Line(1:4) == 'NOSO') then
  NOSO = .true.
  Linenr = Linenr+1
  goto 100
end if
if (Line(1:4) == 'CURD') then
  IFCURD = .true.
  Linenr = Linenr+1
  goto 100
end if
! END BP OPTIONS
! RF SO-NTO
if (line(1:4) == 'SONT') then
  read(LuIn,*,err=997) SONTOSTATES
  call mma_allocate(SONTO,2,SONTOSTATES,Label='SONTO')
  linenr = linenr+1
  read(LuIn,*,err=997) (SONTO(1,I),SONTO(2,I),I=1,SONTOSTATES)
  linenr = linenr+1
  goto 100
end if
if (line(1:4) == 'ARGU') then
  IFARGU = .true.
  Linenr = Linenr+1
  goto 100
end if
! END RF
!-SVC 2007 2008------------------------------
if (Line(1:4) == 'MAGN') then
  IFXCAL = .true.
  read(LuIn,*,err=997) NBSTEP,BSTART,BINCRE,BANGRES
  read(LuIn,*,err=997) NTSTEP,TSTART,TINCRE
  if (BANGRES > Zero) IFMCAL = .true.
  Linenr = Linenr+2
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'XVIN') then
  PRXVR = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'XVES') then
  PRXVE = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'XVSO') then
  PRXVS = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'MEIN') then
  PRMER = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'MEES') then
  PRMEE = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'MESO') then
  PRMES = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
!*****IS 30-09/2007**************************
if (Line(1:4) == 'HOP ') then
  HOP = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'TRAC') then
  TRACK = .true.
  NOHAM = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'STOV') then
  ONLY_OVERLAPS = .true.
  NOHAM = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'DCOU') then
  IfDCpl = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
!*PAM Nov 2011
if (Line(1:4) == 'TRD1') then
  IFTRD1 = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'TRD2') then
  IFTRD1 = .true.
  IFTRD2 = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (Line(1:3) == 'TDM') then
  IFTDM = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
!CEH April 2015
if (Line(1:4) == 'DQVD') then
  DQVD = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (LINE(1:4) == 'ALPH') then
  read(LuIn,*,err=997) ALPHZ
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'BETA') then
  read(LuIn,*,err=997) BETAE
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
!LKS Sep 2015
if (Line(1:4) == 'DIPR') then
  ! Printing threshold for dipole intensities. Current default 1.0e-5
  DIPR = .true.
  read(LuIn,*,err=997) OSTHR_DIPR
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (LINE(1:4) == 'QIPR') then
  ! Printing threshold for quadrupole intensities. Current default 1.0e-5
  QIPR = .true.
  read(LuIn,*,err=997) OSTHR_QIPR
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'QIAL') then
  ! Print all contributions for quadrupole intensities.
  QIALL = .true.
  LINENR = LINENR+1
  goto 100
end if
!--------------------------------------------
if (LINE(1:4) == 'RSPR') then
  ! Printing threshold for rotatory strength. Current default 1.0e-7
  RSPR = .true.
  read(LuIn,*,err=997) RSTHR
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'CD  ') then
  ! Perform regular circular dichroism - velocity and mixed gauge
  DOCD = .true.
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'DYSO') then
  ! Enable Dyson orbital calculations
  DYSO = .true.
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'DYSE') then
  ! Enable Dyson orbital calculations
  DYSEXPORT = .true.
  read(LuIn,*,err=997) DYSEXPSF,DYSEXPSO
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'TDYS') then
  ! Enable 2particle Dyson matrix calculations
  TDYS = .true.
  read(LuIn,*,err=997) OCAN
  do I=1,OCAN
    read(LuIn,'(A)',err=997) OCAA(I)
  end do
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'DCHS') then
  ! Enable computation of DCH intensities
  DCHS = .true.
  read(LuIn,*,err=997) DCHO
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (Line(1:4) == 'TINT') then
  ! Calculate exact isotropically averaged semi-classical intensities
  ! Activate integration of transition intensities
  ! based on the exact non-relativistic Hamiltonian in the weak field
  ! approximation.
  Do_TMOM = .true.
  Linenr = Linenr+1
  goto 100
end if
! ------------------------------------------
if (Line(1:4) == 'TIGR') then
  ! Group exact TINT to reduce computational cost
  ! TMGr_thrs is the tolerance in the relative energy (unitless)
  ! TMGr_thrs only works with SUBS keyword
  read(LuIn,*,err=997) TMGr_thrs
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (LINE(1:4) == 'PRRA') then
  ! Print the raw directions for exact semi-classical intensities
  PRRAW = .true.
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'PRWE') then
  ! Print the weighted directions for exact semi-classical intensities
  PRWEIGHT = .true.
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'TOLE') then
  ! Set tolerance for different gauges - currently 10 percent (0.1)
  ! Defined as Tolerance = ABS(1-O_r/O_p)
  read(LuIn,*,err=997) TOLERANCE
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'SUBS') then
  ! Reduce looping in intensities. Set limit for the inner and outer loop
  REDUCELOOP = .true.
  read(LuIn,*,err=997) LOOPDIVIDE
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (LINE(1:4) == 'NFIN') then
  ! Reduce looping in intensities. Set number of final states
  read(LuIn,*,err=997) LOOPMAX
  LINENR = LINENR+1
  goto 100
end if
! ------------------------------------------
if (Line(1:4) == 'IIOR') then
  ! Set the order of the Lebedev polynomials used for the numerical
  ! isotropic integration. Current default 5.
  read(LuIn,*,err=997) L_Eff
  ! make sure it's an odd number
  ! and find the smallest grid that supports that degree
  if (mod(L_Eff,2) == 0) L_Eff = L_Eff+1
  do i=(L_Eff-1)/2,rule_max
    if (available_table(i) == 1) exit
    L_Eff = L_Eff+2
  end do
  if (L_Eff > rule_max) then
    call WarningMessage(2,'L_Eff too large, grid not supported')
    call abend()
  end if
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'DIRE') then
  ! Set a specific direction of the incident light when computing
  ! the transition intensities in the use of
  ! the vector field (A) in the non-relativistic Hamiltonian.
  Do_SK = .true.
  read(LuIn,*,err=997) nk_Vector
  Linenr = Linenr+1
  call mma_allocate(k_Vector,3,nk_Vector,label='k-Vector')
  do j=1,nk_Vector
    read(LuIn,*,err=997) (k_Vector(i,j),i=1,3)
    Linenr = Linenr+1
    tmp = k_Vector(1,j)**2+k_Vector(2,j)**2+k_Vector(3,j)**2
    tmp = One/sqrt(tmp)
    k_Vector(1,j) = k_Vector(1,j)*tmp
    k_Vector(2,j) = k_Vector(2,j)*tmp
    k_Vector(3,j) = k_Vector(3,j)*tmp
  end do
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'POLA') then
  Do_Pol = .true.
  Linenr = Linenr+1
  read(LuIn,*,err=997) (e_Vector(i),i=1,3)
  goto 100
end if
!--------------------------------------------
! VKochetov 2021 enable saving more data to hdf5
if (Line(1:4) == 'RHOD') then
  rhodyn = .true.
  Linenr = Linenr+1
  goto 100
end if
!--------------------------------------------
#ifdef _DMRG_
if (Line(1:4) == 'QDSC') then
  QDPT2SC = .true.
  goto 100
end if
!--------------------------------------------
if (Line(1:4) == 'QDPC') then
  QDPT2SC = .false.
  goto 100
end if
#endif
! These errors could eventually be removed
!--------------------------------------------
if (Line(1:4) == 'TMOS') then
  write(u6,*) ' Keyword TMOS is deprecated, use TINT instead.'
  goto 999
end if
!--------------------------------------------
if (Line(1:4) == 'KVEC') then
  write(u6,*) ' Keyword KVEC is deprecated, use DIRE instead.'
  goto 999
end if
!--------------------------------------------
if (Line(1:4) == 'L-EF') then
  write(u6,*) ' Keyword L-EF is deprecated, use IIOR instead.'
  goto 999
end if
!--------------------------------------------
if (Line(1:4) == 'REDL') then
  write(u6,*) ' Keyword REDL is deprecated, use SUBS instead.'
  goto 999
end if
!--------------------------------------------
if (Line(1:4) == 'TMGR') then
  write(u6,*) ' Keyword TMGR is deprecated, use TIGR instead.'
  goto 999
end if
!--------------------------------------------
write(u6,*) ' The following input line was not understood:'
write(u6,'(A)') LINE
goto 999

997 continue
call WarningMessage(2,'Error reading standard input.')
write(u6,*) ' RASSI input near line nr.',LINENR+1
goto 999

998 continue
call WarningMessage(2,'I/O error.')
write(u6,*) ' READIN: Unexpected end of input file.'

999 continue
call ABEND()

200 continue
!nf
if (IfDCpl .and. (.not. IfHam)) then
  call WarningMessage(1,'Input request was ignored.')
  write(u6,*) ' Cannot compute the approximate derivative coupling terms without the energies.'
  write(u6,*) ' Ignore them.'
  IfDCpl = .false.
end if
!nf
if (Do_Pol .and. (.not. Do_SK)) then
  call WarningMessage(1,'Input request was ignored.')
  write(u6,*) ' Polarization direction can only be used with'
  write(u6,*) ' specific k-vector directions.'
  Do_Pol = .false.
end if
! Prints warning if rot. str. threshold is defined without any calculations
if (RSPR) then
  if ((.not. DOCD) .and. (.not. Do_TMOM)) then
    call WarningMessage(1,'Input request was ignored.')
    write(u6,*) 'Warning: Rotatory strength threshold specified (RSPR) without calculating rotatory strength'
  end if
end if
! Determine file names, if undefined.
if (JBNAME(1) == 'UNDEFINE') then
  ! The first (perhaps only) jobiph file is named 'JOB001', or maybe 'JOBIPH'
  ! when no name has been issued by the user:
  NFLS = 0
  TRYNAME = 'JOB001'
  call F_INQUIRE(TRYNAME,LEXISTS)
  if (LEXISTS) then
    NFLS = 1
    JBNAME(NFLS) = TRYNAME
  else
    TRYNAME = 'JOBIPH'
    call F_INQUIRE(TRYNAME,LEXISTS)
    if (LEXISTS) then
      NFLS = 1
      JBNAME(NFLS) = TRYNAME
    else
      call WarningMessage(1,'RASSI lacks JobIph files.')
      write(u6,*) ' RASSI fails: No jobiph files found.'
      call ABEND()
    end if
  end if
  ! Subsequent (if any) jobfiles can be named according to old
  ! or new naming convention.
  ! Using new standard scheme for default jobiph names?
  do I=1,MXJOB-1
    if ((NJOB > 0) .and. (I > NJOB)) goto 211
    write(TRYNAME,'(A6,I2.2)') 'JOBIPH',I
    call F_INQUIRE(TRYNAME,LEXISTS)
    if (LEXISTS) then
      NFLS = NFLS+1
      JBNAME(NFLS) = TRYNAME
    else
      goto 211
    end if
  end do
  call WarningMessage(1,'RASSI fails to identify JobIph files.')
  write(u6,*) ' Too many jobiph files in this directory.'
  call ABEND()
211 continue
  if (NFLS == 1) then
    ! We may be using old standard scheme for default jobiph names?
    do I=1,MXJOB-1
      if ((NJOB > 0) .and. (I > NJOB)) goto 212
      write(TRYNAME,'(A3,I3.3)') 'JOB',I+1
      call F_INQUIRE(TRYNAME,LEXISTS)
      if (LEXISTS) then
        NFLS = NFLS+1
        JBNAME(NFLS) = TRYNAME
      else
        goto 212
      end if
    end do
    call WarningMessage(1,'RASSI fails to identify JobIph files.')
    write(u6,*) ' Too many jobiph files in this directory.'
    call ABEND()
212 continue
    ! Then we are definitely using the old default file name convention.
    !if (NFLS > 1) then
    !end if
  end if
  if (NJOB > 0) then
    ! Input has been given for NJOB, etc., and will be used.
    if (NFLS < NJOB) then
      call WarningMessage(1,'RASSI found too few JobIph files.')
      call ABEND()
    end if
  else
    ! Use defaults.
    NJOB = NFLS
  end if
end if

call XFLUSH(u6)

call Close_LuSpool(LuIn)

return

end subroutine READIN_RASSI
