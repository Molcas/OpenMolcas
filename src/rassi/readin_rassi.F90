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

#include "macros.fh"

subroutine READIN_RASSI()

use Cholesky, only: timings
use Cntrl, only: ALGO, Nscreen, dmpk, QDPT2SC, QDPT2EV, SECOND_TIME, DOGSOR, PRSXY, PRORB, PRTRA, PRCI, BINA, NATO, NBINA, NRNATO, &
                 RFPERT, IFTRD1, NSOPR, NPROP, PRDIPVEC, TDIPMIN, NJOB, CIH5, CIThr, IFHAM, IFSO, IFNTO, SOThr_Prt, nSOThr_Prt, &
                 nState, IfHEXT, IfHEff, IfHCOM, IFEJOB, IfHDia, IfShft, ToFile, IfJ2, IfJZ, IFGCAL, EPraThr, IFACALSD, IFACALFC, &
                 IFACALSDON, IFACALPSO, IFATCALSA, IFGTSHSA, MULTIP, IFVANVLECK, TMINS, TMAXS, NTS, IFSONCINI, TMINP, TMAXP, NTP, &
                 NOSO, IFCURD, IFARGU, IFXCAL, NBSTEP, BSTART, BINCRE, BANGRES, NTSTEP, TSTART, TINCRE, IFMCAL, PRXVR, PRXVE, &
                 PRXVS, PRMER, PRMEE, PRMES, HOP, TRACK, NOHAM, ONLY_OVERLAPS, IFDCPL, IFTRD2, IFTDM, DQVD, ALPHZ, BETAE, DIPR, &
                 OSTHR_DIPR, QIPR, OSTHR_QIPR, QIALL, RSPR, RSThr, DOCD, DYSO, DYSEXPORT, DYSEXPSO, TDYS, OCAN, DCHS, DCHO, &
                 DO_TMOM, TMGR_Thrs, PRRAW, PRWEIGHT, TOLERANCE, REDUCELOOP, LOOPDIVIDE, LOOPMAX, l_Eff, Do_SK, Do_Pol, RHODYN, &
                 MXJOB, JBNAME, SOPRNM, PNAME, PRDIPCOM, EPrThr, LPRPR, lHami, IfACAL, IFACALFCON, IFACALFCSDON, IFGTCALSA, &
                 DYSEXPSF, ISTAT, MXPROP, NSTAT, IBINA, ISOCMP, ICOMP, OCAA, SONTO, SONTOSTATES, SONAT, SONATNSTATE, SODIAG, &
                 SODIAGNSTATE
use Fock_util_global, only: Deco, Estimate, PseudoChoMOs, Update
use frenkel_global_vars, only: DoCoul, doexch, DoExcitonics, excl, iTyp, labB, nestla, nestlb, valst
use kVectors, only: e_Vector, k_Vector, nk_Vector
use Lebedev_quadrature, only: available_table, rule_max
use rassi_data, only: CHFRACMEM
use rassi_global_arrays, only: ESHFT, HAM, HDIAG, JBNUM, LROOT
use spool, only: Close_LuSpool, Spoolinp
#ifdef _DMRG_
use rasscf_global, only: doDMRG
#endif
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: corest, I, IJOB, ISTATE, istatus, J, JSTATE, LINENR, LuIn, nesta, nestb, NFLS
real(kind=wp) :: tmp
logical(kind=iwp) :: lExists
character(len=80) ::LINE
character(len=8) :: TRYNAME
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
LINENR = 0
do
  read(LuIn,'(A72)',iostat=istatus) LINE
  call LineCheck(istatus)
  call NORMAL(LINE)
  if (LINE(1:7) == input_id) exit
end do
LINENR = 0
do
  read(LuIn,'(A72)',iostat=istatus) LINE
  call LineCheck(istatus)
  call NORMAL(LINE)
  if (LINE(1:1) == '*') cycle
  if (LINE == '') cycle
  select case (LINE(1:4))

    case ('END ')
      exit

    case ('TEST')
      PRSXY = .true.
      PRORB = .true.
      PRTRA = .true.
      PRCI = .true.

    case ('SECO')
      SECOND_TIME = .true.

    case ('GSOR')
      DoGSOR = .true.

    case ('BINA')
      BINA = .true.
      NATO = .true.
      read(LuIn,*,iostat=istatus) NBINA
      call LineCheck(istatus)
      read(LuIn,*,iostat=istatus) (IBINA(1,I),IBINA(2,I),I=1,NBINA)
      call LineCheck(istatus)

    case ('NATO')
      NATO = .true.
      read(LuIn,*,iostat=istatus) NRNATO
      call LineCheck(istatus)

    case ('RFPE')
      RFPERT = .true.

    case ('CHOL')
      ! --- Cholesky with default settings
      call Cho_rassi_rdInp(.true.,LuIn)

    case ('CHOI')
      ! --- Cholesky with customized settings
      call Cho_rassi_rdInp(.false.,LuIn)

    case ('EXAL')
      EXCL = .true.
      NESTA = 0
      read(LuIn,*,iostat=istatus) NESTA
      call LineCheck(istatus)
      call mma_allocate(NESTLA,NESTA)
      read(LuIn,*,iostat=istatus) (NESTLA(I),I=1,NESTA)
      call LineCheck(istatus)
      !write(u6,*) 'list of excit. initial states'
      !write(u6,*) 'mon A nr states', NESTA
      !write(u6,'(I2)') (NESTLA(I),I=1,NESTA)

    case ('EXBL')
      EXCL = .true.
      NESTB = 0
      read(LuIn,*,iostat=istatus) NESTB
      call LineCheck(istatus)
      call mma_allocate(NESTLB,NESTB)
      read(LuIn,*,iostat=istatus) (NESTLB(I),I=1,NESTB)
      call LineCheck(istatus)
      !write(u6,*) 'list of excit. initial states'
      !write(u6,*) 'mon B nr states', NESTB
      !write(u6,'(I2)') (NESTLB(I),I=1,NESTB)

    case ('EXCI')
      DoExcitonics = .true.

    case ('KCOU')
      DoExch = .true.

    case ('MONA')
      DoCoul = .true.
      if (iTyp == 2) write(u6,*) ' Warning: switching to monomer-A.'
      iTyp = 1
      if (.not. IFTRD1) then
        IFTRD1 = .true.
        write(u6,*) ' TRD1 activated by MONA.'
      end if

    case ('MONB')
      labB = .true.
      DoCoul = .true.
      if (iTyp == 1) write(u6,*) ' Warning: switching to monomer-B.'
      iTyp = 2
      if (.not. IFTRD1) then
        IFTRD1 = .true.
        write(u6,*) ' TRD1 activated by MONB.'
      end if

    case ('RIXS')
      read(LuIn,*,iostat=istatus) valst,corest
      call LineCheck(istatus)
      unused_var(corest)

    case ('SOPR')
      read(LuIn,*,iostat=istatus) NSOPR,(SOPRNM(I),ISOCMP(I),I=1,min(MXPROP,NSOPR))
      call LineCheck(istatus)
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

    case ('PROP')
      read(LuIn,*,iostat=istatus) NPROP,(PNAME(I),ICOMP(I),I=1,min(MXPROP,NPROP))
      call LineCheck(istatus)
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

    case ('OVER')
      PRSXY = .true.

    case ('TRDI')
      ! Print transition dipole vectors
      PRDIPVEC = .true.

    case ('TDMN')
      ! Print transition dipole vectors
      PRDIPVEC = .true.
      read(LuIn,*,iostat=istatus) TDIPMIN
      call LineCheck(istatus)

    case ('TRDC')
      ! Print COMPLEX transition dipole vectors
      PRDIPCOM = .true.

    case ('ORBI')
      PRORB = .true.

    case ('CIPR')
      PRCI = .true.

    case ('CIH5')
      if (NJOB <= 2) then
        CIH5 = .true.
      else
        call WarningMessage(2,'CIH5 allows no more than 2 JOBIPHs')
        call abend()
      end if

    case ('THRS')
      read(LuIn,*,iostat=istatus) CITHR
      call LineCheck(istatus)

    case ('ONEL')
      IFHAM = .false.

    case ('ONEE')
      IFHAM = .false.

    case ('SPIN')
      IFSO = .true.

    case ('NTOC')
      IFNTO = .true.

    case ('SOCO')
      ! PAM07 Added: Keyword for printing spin-orbit coupling matrix elements
      ! A threshold in reciprocal cm is entered.
      read(LuIn,*,iostat=istatus) SOTHR_PRT
      call LineCheck(istatus)
      if (SOTHR_PRT < Zero) SOTHR_PRT = Zero
      NSOTHR_PRT = 10000

    case ('NROF','NR O')
      NSTATE = 0
      read(LuIn,*,iostat=istatus) NJOB,TRYNAME
      call LineCheck(istatus)
      call UpCase(TRYNAME)
      if (TRYNAME /= 'ALL') then
        backspace(LuIn)
        LINENR = LINENR-1
        read(LuIn,*,iostat=istatus) NJOB,(NSTAT(I),I=1,NJOB)
        call LineCheck(istatus)
        do IJOB=1,NJOB
          NSTATE = NSTATE+NSTAT(IJOB)
        end do
        call mma_allocate(JBNUM,nState,Label='JBNUM')
        call mma_allocate(LROOT,nState,Label='LROOT')
        NSTATE = 0
        do IJOB=1,NJOB
          ISTAT(IJOB) = NSTATE+1
          read(LuIn,*,iostat=istatus) (LROOT(NSTATE+J),J=1,NSTAT(IJOB))
          call LineCheck(istatus)
          do ISTATE=NSTATE+1,NSTATE+NSTAT(IJOB)
            JBNUM(ISTATE) = IJOB
          end do
          NSTATE = NSTATE+NSTAT(IJOB)
        end do
      end if

    case ('IPHN')
      !read(LuIn,'(9(A7,1X))',iostat=istatus) (JBNAME(I),I=1,NJOB)
      !call LineCheck(istatus)
      do I=1,NJOB
        read(LuIn,*,iostat=istatus) JBNAME(I)
        call LineCheck(istatus)
        call UpCase(JBNAME(I))
      end do

    case ('FILE')
      read(LuIn,*,iostat=istatus) NJOB
      call LineCheck(istatus)
      do I=1,NJOB
        read(LuIn,*,iostat=istatus) LINE
        call LineCheck(istatus)
        call FILEORB(LINE,JBNAME(I))
      end do

    case ('HEXT')
      if (nState == 0) then
        call WarningMessage(2,'HEXT needs an explicit NROFJOB')
        call abend()
      end if
      IFHEXT = .true.
      call mma_allocate(HAM,nState,nState,Label='HAM')
      read(LuIn,*,iostat=istatus) ((HAM(ISTATE,JSTATE),JSTATE=1,ISTATE),ISTATE=1,NSTATE)
      call LineCheck(istatus)
      do ISTATE=1,NSTATE-1
        do JSTATE=ISTATE+1,NSTATE
          HAM(ISTATE,JSTATE) = HAM(JSTATE,ISTATE)
        end do
      end do

    case ('HEFF')
      IFHEFF = .true.

    case ('HCOM')
      IFHCOM = .true.

    case ('EJOB')
      IFEJOB = .true.

    case ('HDIA')
      if (nState == 0) then
        call WarningMessage(2,'HDIA needs an explicit NROFJOB')
        call abend()
      end if
      IFHDIA = .true.
      call mma_allocate(HDIAG,nState,Label='nState')
      read(LuIn,*,iostat=istatus) (HDIAG(ISTATE),ISTATE=1,NSTATE)
      call LineCheck(istatus)

    case ('SHIF')
      if (nState == 0) then
        call WarningMessage(2,'SHIFT needs an explicit NROFJOB')
        call abend()
      end if
      IFSHFT = .true.
      call mma_allocate(ESHFT,nState,Label='ESHFT')
      read(LuIn,*,iostat=istatus) (ESHFT(ISTATE),ISTATE=1,NSTATE)
      call LineCheck(istatus)

    case ('TOFI')
      ToFile = .true.

    case ('J-VA')
      IFJ2 = 1

    case ('OMEG')
      IFJZ = 1

    case ('EPRG')
      !-SVC 2007-----------------------------------
      IFGCAL = .true.
      read(LuIn,*,iostat=istatus) EPRTHR
      call LineCheck(istatus)
      if (EPRTHR < Zero) EPRTHR = Zero

    case ('PRPR')
      ! tjd- BMII: Print out spin-orbit properties to files
      write(u6,*) 'SPIN-ORBIT PROPERTY PRINT ON'
      LPRPR = .true.

    case ('HAMI')
      ! tjd- Yoni: Force an identity SO hamiltonian
      write(u6,*) 'Identity Hamiltonian turned on'
      LHAMI = .true.

    case ('EPRA')
      ! BP - Hyperfine calculations
      !write(u6,*)"EPRA read"
      IFACAL = .true.
      read(LuIn,*,iostat=istatus) EPRATHR
      call LineCheck(istatus)
      if (EPRATHR < Zero) EPRATHR = Zero

    case ('AFCO')
      !write(u6,*) 'AFCO read'
      IFACALSD = .false.

    case ('ASDO')
      !write(u6,*) 'ASDO read'
      IFACALFC = .false.

    case ('AFCC')
      ! Kamal Sharkas beg - PSO Hyperfine calculations
      IFACALFCON = .true.

    case ('ASDC')
      IFACALSDON = .true.

    case ('FCSD')
      IFACALFCSDON = .true.

    case ('APSO')
      IFACALPSO = .true.

    case ('GTSA')
      IFGTCALSA = .true.

    case ('ATSA')
      IFATCALSA = .true.

    case ('SHMP')
      IFGTSHSA = .true.
      read(LuIn,*,iostat=istatus) MULTIP
      call LineCheck(istatus)

    case ('VSUS')
      IFVANVLECK = .true.
      read(LuIn,*,iostat=istatus) TMINS,TMAXS,NTS
      call LineCheck(istatus)

    case ('NMRT')
      IFSONCINI = .true.
      read(LuIn,*,iostat=istatus) TMINP,TMAXP,NTP
      call LineCheck(istatus)
      ! Kamal Sharkas end - PSO Hyperfine calculations

    case ('SONO')
      ! BP Natural orbitals options
      read(LuIn,*,iostat=istatus) SONATNSTATE
      call LineCheck(istatus)
      call mma_allocate(SONAT,SONATNSTATE,Label='SONAT')
      read(LuIn,*,iostat=istatus) (SONAT(I),I=1,SONATNSTATE)
      call LineCheck(istatus)

    case ('SODI')
      read(LuIn,*,iostat=istatus) SODIAGNSTATE
      call LineCheck(istatus)
      call mma_allocate(SODIAG,SODIAGNSTATE,Label='SODIAG')
      read(LuIn,*,iostat=istatus) (SODIAG(I),I=1,SODIAGNSTATE)
      call LineCheck(istatus)

    case ('NOSO')
      NOSO = .true.

    case ('CURD')
      IFCURD = .true.
      ! END BP OPTIONS

    case ('SONT')
      ! RF SO-NTO
      read(LuIn,*,iostat=istatus) SONTOSTATES
      call LineCheck(istatus)
      call mma_allocate(SONTO,2,SONTOSTATES,Label='SONTO')
      read(LuIn,*,iostat=istatus) (SONTO(1,I),SONTO(2,I),I=1,SONTOSTATES)
      call LineCheck(istatus)

    case ('ARGU')
      IFARGU = .true.
      ! END RF

    case ('MAGN')
      IFXCAL = .true.
      read(LuIn,*,iostat=istatus) NBSTEP,BSTART,BINCRE,BANGRES
      call LineCheck(istatus)
      read(LuIn,*,iostat=istatus) NTSTEP,TSTART,TINCRE
      call LineCheck(istatus)
      if (BANGRES > Zero) IFMCAL = .true.

    case ('XVIN')
      PRXVR = .true.

    case ('XVES')
      PRXVE = .true.

    case ('XVSO')
      PRXVS = .true.

    case ('MEIN')
      PRMER = .true.

    case ('MEES')
      PRMEE = .true.

    case ('MESO')
      PRMES = .true.

    case ('HOP ')
      HOP = .true.

    case ('TRAC')
      TRACK = .true.
      NOHAM = .true.

    case ('STOV')
      ONLY_OVERLAPS = .true.
      NOHAM = .true.

    case ('DCOU')
      IfDCpl = .true.

    case ('TRD1')
      !*PAM Nov 2011
      IFTRD1 = .true.

    case ('TRD2')
      IFTRD1 = .true.
      IFTRD2 = .true.

    case ('TDM ')
      IFTDM = .true.

    case ('DQVD')
      !CEH April 2015
      DQVD = .true.

    case ('ALPH')
      read(LuIn,*,iostat=istatus) ALPHZ
      call LineCheck(istatus)

    case ('BETA')
      read(LuIn,*,iostat=istatus) BETAE
      call LineCheck(istatus)

    case ('DIPR')
      !LKS Sep 2015
      ! Printing threshold for dipole intensities. Current default 1.0e-5
      DIPR = .true.
      read(LuIn,*,iostat=istatus) OSTHR_DIPR
      call LineCheck(istatus)

    case ('QIPR')
      ! Printing threshold for quadrupole intensities. Current default 1.0e-5
      QIPR = .true.
      read(LuIn,*,iostat=istatus) OSTHR_QIPR
      call LineCheck(istatus)

    case ('QIAL')
      ! Print all contributions for quadrupole intensities.
      QIALL = .true.

    case ('RSPR')
      ! Printing threshold for rotatory strength. Current default 1.0e-7
      RSPR = .true.
      read(LuIn,*,iostat=istatus) RSTHR
      call LineCheck(istatus)

    case ('CD  ')
      ! Perform regular circular dichroism - velocity and mixed gauge
      DOCD = .true.

    case ('DYSO')
      ! Enable Dyson orbital calculations
      DYSO = .true.

    case ('DYSE')
      ! Enable Dyson orbital calculations
      DYSEXPORT = .true.
      read(LuIn,*,iostat=istatus) DYSEXPSF,DYSEXPSO
      call LineCheck(istatus)

    case ('TDYS')
      ! Enable 2particle Dyson matrix calculations
      TDYS = .true.
      read(LuIn,*,iostat=istatus) OCAN
      call LineCheck(istatus)
      do I=1,OCAN
        read(LuIn,'(A)',iostat=istatus) OCAA(I)
        call LineCheck(istatus)
      end do

    case ('DCHS')
      ! Enable computation of DCH intensities
      DCHS = .true.
      read(LuIn,*,iostat=istatus) DCHO
      call LineCheck(istatus)

    case ('TINT')
      ! Calculate exact isotropically averaged semi-classical intensities
      ! Activate integration of transition intensities
      ! based on the exact non-relativistic Hamiltonian in the weak field
      ! approximation.
      Do_TMOM = .true.

    case ('TIGR')
      ! Group exact TINT to reduce computational cost
      ! TMGr_thrs is the tolerance in the relative energy (unitless)
      ! TMGr_thrs only works with SUBS keyword
      read(LuIn,*,iostat=istatus) TMGr_thrs
      call LineCheck(istatus)

    case ('PRRA')
      ! Print the raw directions for exact semi-classical intensities
      PRRAW = .true.

    case ('PRWE')
      ! Print the weighted directions for exact semi-classical intensities
      PRWEIGHT = .true.

    case ('TOLE')
      ! Set tolerance for different gauges - currently 10 percent (0.1)
      ! Defined as Tolerance = ABS(1-O_r/O_p)
      read(LuIn,*,iostat=istatus) TOLERANCE
      call LineCheck(istatus)

    case ('SUBS')
      ! Reduce looping in intensities. Set limit for the inner and outer loop
      REDUCELOOP = .true.
      read(LuIn,*,iostat=istatus) LOOPDIVIDE
      call LineCheck(istatus)

    case ('NFIN')
      ! Reduce looping in intensities. Set number of final states
      read(LuIn,*,iostat=istatus) LOOPMAX
      call LineCheck(istatus)

     case ('IIOR')
      ! Set the order of the Lebedev polynomials used for the numerical
      ! isotropic integration. Current default 5.
      read(LuIn,*,iostat=istatus) L_Eff
      call LineCheck(istatus)
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

    case ('DIRE')
      ! Set a specific direction of the incident light when computing
      ! the transition intensities in the use of
      ! the vector field (A) in the non-relativistic Hamiltonian.
      Do_SK = .true.
      read(LuIn,*,iostat=istatus) nk_Vector
      call LineCheck(istatus)
      call mma_allocate(k_Vector,3,nk_Vector,label='k-Vector')
      do j=1,nk_Vector
        read(LuIn,*,iostat=istatus) (k_Vector(i,j),i=1,3)
        call LineCheck(istatus)
        tmp = k_Vector(1,j)**2+k_Vector(2,j)**2+k_Vector(3,j)**2
        tmp = One/sqrt(tmp)
        k_Vector(1,j) = k_Vector(1,j)*tmp
        k_Vector(2,j) = k_Vector(2,j)*tmp
        k_Vector(3,j) = k_Vector(3,j)*tmp
      end do

    case ('POLA')
      Do_Pol = .true.
      read(LuIn,*,iostat=istatus) (e_Vector(i),i=1,3)
      call LineCheck(istatus)

    case ('RHOD')
      ! VKochetov 2021 enable saving more data to hdf5
      rhodyn = .true.

#   ifdef _DMRG_
    case ('QDSC')
      QDPT2SC = .true.

    case ('QDPC')
      QDPT2SC = .false.
#   endif

    case ('EXTR')
      call WarningMessage(2,' Keyword EXTR is obsolete.')

    case ('TMOS')
      call WarningMessage(2,' Keyword TMOS is deprecated, use TINT instead.')

    case ('KVEC')
      call WarningMessage(2,' Keyword KVEC is deprecated, use DIRE instead.')

    case ('L-EF')
      call WarningMessage(2,' Keyword L-EF is deprecated, use IIOR instead.')

    case ('REDL')
      call WarningMessage(2,' Keyword REDL is deprecated, use SUBS instead.')

    case ('TMGR')
      call WarningMessage(2,' Keyword TMGR is deprecated, use TIGR instead.')

    case default
      write(u6,*) ' The following input line was not understood:'
      write(u6,'(A)') LINE
      call abend()

  end select

end do

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
    if ((NJOB > 0) .and. (I > NJOB)) exit
    write(TRYNAME,'(A6,I2.2)') 'JOBIPH',I
    call F_INQUIRE(TRYNAME,LEXISTS)
    if (LEXISTS) then
      NFLS = NFLS+1
      JBNAME(NFLS) = TRYNAME
    else
      exit
    end if
  end do
  if (I > MXJOB-1) then
    call WarningMessage(1,'RASSI fails to identify JobIph files.')
    write(u6,*) ' Too many jobiph files in this directory.'
    call ABEND()
  end if
  if (NFLS == 1) then
    ! We may be using old standard scheme for default jobiph names?
    do I=1,MXJOB-1
      if ((NJOB > 0) .and. (I > NJOB)) exit
      write(TRYNAME,'(A3,I3.3)') 'JOB',I+1
      call F_INQUIRE(TRYNAME,LEXISTS)
      if (LEXISTS) then
        NFLS = NFLS+1
        JBNAME(NFLS) = TRYNAME
      else
        exit
      end if
    end do
    if (I > MXJOB-1) then
      call WarningMessage(1,'RASSI fails to identify JobIph files.')
      write(u6,*) ' Too many jobiph files in this directory.'
      call ABEND()
    end if
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

contains

subroutine LineCheck(code)

  integer(kind=iwp), intent(in) :: code

  LineNr = LineNr+1
  select case (code)
    case (:-1)
      call WarningMessage(2,'I/O error.')
      write(u6,*) ' READIN: Unexpected end of input file.'
      call ABEND()
    case (0)
    case (1:)
      call WarningMessage(2,'Error reading standard input.')
      write(u6,*) ' RASSI input near line nr.',LINENR
      call ABEND()
  end select

end subroutine LineCheck

end subroutine READIN_RASSI
