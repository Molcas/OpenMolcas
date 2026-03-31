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

subroutine INIT_RASSI()

use symmetry_Info, only: symmetry_info_get
use rassi_aux, only: ipglob
#ifndef _DMRG_
use rasscf_global, only: doDMRG
#endif
use Cntrl, only: BINA, CIH5, CITHR, DCHO, DCHS, DIPR, Do_Pol, Do_SK, DO_TMOM, DoCD, DYSEXPORT, DYSO, EPRAThr, EPRThr, FnEig, &
                 FnTOM, FORCE_NON_AO_TDM, HAVE_DIAG, HAVE_HEFF, HOP, IFACAL, IFACALFC, IFACALSD, IFArgU, IFATCALSA, IFCURD, &
                 IFDCPL, IFEJOB, IFGCAL, IFGTCALSA, IFGTSHSA, IFHAM, IFHCOM, IFHDIA, IFHEFF, IFHEXT, IFMCAL, IFNTO, IFSHFT, IFSO, &
                 IFTDM, IFTRD1, IFTRD2, IFXCAL, JBNAME, L_Eff, LHAMI, LOOPDIVIDE, LOOPMAX, LPRPR, LuEig, LuExc, LuIph, LuMck, &
                 LuOne, LuOrd, LuTOM, MINAME, MXPROP, NATO, NBINA, NJOB, NOHAM, NOSO, NPROP, NrNATO, NSOPR, NSOThr_Prt, NSTATE, &
                 OCAN, ONLY_OVERLAPS, OSThr_DipR, OSThr_QIPR, PNAME, PRCI, PRDIPVEC, PRMEE, PRMER, PRMES, PRORB, PRRAW, PRSXY, &
                 PRTRA, PRWEIGHT, PRXVE, PRXVR, PRXVS, PTYPE, QIALL, QIPR, REDUCELOOP, RFPert, RSPR, RSThr, SODIAGNSTATE, &
                 SONATNSTATE, SONTOSTATES, SOPRNM, SOPRTP, SOThr_Prt, TDIPMIN, TDYS, TMGR_Thrs, ToFile, TOLERANCE, TRACK
use rassi_data, only: WFTYPE
use hfc_logical, only: MAG_X2C
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, IPROP
logical(kind=iwp) :: DoCholesky, FoundTwoEls
character(len=256) :: STRING

! Initialise doDMRG if compiled without QCMaquis
#ifndef _DMRG_
DoDMRG = .false.
#endif

call symmetry_info_get()

! UNIT NUMBERS AND NAMES
LUONE = 2
LUORD = 30
LUIPH = 15
LUEXC = 22
LUMCK = 33
LuToM = 26
FnToM = 'TOFILE'
LuEig = 27
FnEig = 'EIGV'
JBNAME(:) = 'UNDEFINE'
do I=1,size(JBNAME)
  write(MINAME(I),'(''MCK'',I3.3)') I
end do

if (IPGLOB > 3) write(u6,*) ' OPENING ','ANNI'
call DANAME(LUEXC,'ANNI')

! NR OF JOBIPHS AND STATES:
NJOB = 0
NSTATE = 0
if (IPGLOB > 3) then
  write(u6,*) ' INITIAL DEFAULT VALUES:'
  write(u6,'(1X,A,I4)') '  NJOB:',NJOB
  write(u6,'(1X,A,I4)') 'NSTATE:',NSTATE
end if

! NR OF OPERATORS FOR WHICH MATRIX ELEMENTS ARE TO BE CALCULATED:
NPROP = 0

! OPERATORS FOR WHICH MATRIX ELEMENTS OVER SPIN-ORBIT EIGENSTATES
! ARE TO BE COMPUTED.
NSOPR = 0

! DEFAULT THRESHOLD FOR PRINTING CI COEFFICIENTS:
CITHR = 0.05_wp

! DEFAULT THRESHOLD FOR PRINTING TRANSITION DIPOLE VECTORS
TDIPMIN = 1.0e-4_wp

! DEFAULT THRESHOLD AND MAX NUMBER OF SO-HAMILTONIAN
! MATRIX ELEMENTS TO PRINT:
NSOTHR_PRT = 0
SOTHR_PRT = -One

! SET LABELS TO UNDEFINED
do IPROP=1,MXPROP
  PNAME(IPROP) = 'UNDEF.'
  PTYPE(IPROP) = 'UNDEF.'
  SOPRNM(IPROP) = 'UNDEF.'
  SOPRTP(IPROP) = 'UNDEF.'
end do

! DEFAULT FLAGS:
PRSXY = .false.
PRDIPVEC = .false.
PRORB = .false.
PRTRA = .false.
PRCI = .false.
CIH5 = .false.
IFHAM = .false.
IFEJOB = .false.
IFSHFT = .false.
IFHDIA = .false.
IFHEXT = .false.
IFHEFF = .false.
IFHCOM = .false.
HAVE_HEFF = .false.
HAVE_DIAG = .false.
NOHAM = .false.
IFSO = .false.
IFNTO = .false.
NATO = .false.
BINA = .false.
IFTRD1 = .false.
IFTRD2 = .false.
IFTDM = .false.
RFPERT = .false.
ToFile = .false.
PRXVR = .false.
PRXVE = .false.
PRXVS = .false.
PRMER = .false.
PRMEE = .false.
PRMES = .false.
IFGCAL = .false.
EPRTHR = Zero
EPRATHR = Zero
IFXCAL = .false.
IFMCAL = .false.
HOP = .false.
TRACK = .false.
ONLY_OVERLAPS = .false.
! Intensities
DIPR = .false.
OSTHR_DIPR = Zero
QIPR = .false.
OSTHR_QIPR = Zero
QIALL = .false.
DYSO = .false.
DYSEXPORT = .false.
TDYS = .false.
OCAN = 1
DCHS = .false.
DCHO = 1
! Exact operator
Do_TMOM = .false.
PRRAW = .false.
PRWEIGHT = .false.
TOLERANCE = 0.1_wp
REDUCELOOP = .false.
LOOPDIVIDE = 0
LOOPMAX = -1
TMGr_thrs = -One
Do_SK = .false.
Do_Pol = .false.
L_Eff = 5
! CD - velocity and mixed gauge
DOCD = .false.
RSTHR = Zero
RSPR = .false.
! Force that TDMs are not stored in the AO basis.
Force_NON_AO_TDM = .false.
call GETENVF('MOLCAS_FORCE_NON_AO_TDM',STRING)
if (STRING == 'ON') Force_NON_AO_TDM = .true.
!nf
IfDCpl = .false.
!nf

! tjd- BMII: Print out spin-orbit properties to files
LPRPR = .false.
LHAMI = .false.
! Feng: test control
MAG_X2C = .false.

! K. Sharkas  BEG
IFATCALSA = .false.
IFGTCALSA = .false.
IFGTSHSA = .false.
! K. Sharkas  END

! BP - Hyperfine tensor and SONATORB initialization
! RF - SO-NTO initialization
IFACAL = .false.
IFACALFC = .true.
IFACALSD = .true.

NOSO = .false.
SONATNSTATE = 0
SODIAGNSTATE = 0

SONTOSTATES = 0

IFCURD = .false.
IFARGU = .false.

! Nr of states for which natural orbitals will be computed:
NRNATO = 0
! Nr of state pairs for computing bi-natural orbitals:
NBINA = 0

! Check if two-electron integrals are available:
call f_Inquire('ORDINT',FoundTwoEls)
call DecideOnCholesky(DoCholesky)
if (FoundTwoEls .or. DoCholesky) IFHAM = .true.

if (IPGLOB >= 4) then
  write(u6,*) 'Initial default flags are:'
  write(u6,*) '     PRSXY :',PRSXY
  write(u6,*) '     PRORB :',PRORB
  write(u6,*) '     PRTRA :',PRTRA
  write(u6,*) '     PRCI  :',PRCI
  write(u6,*) '     IFHAM :',IFHAM
  write(u6,*) '     IFHEXT:',IFHEXT
  write(u6,*) '     IFHEFF:',IFHEFF
  write(u6,*) '     IFEJOB:',IFEJOB
  write(u6,*) '     IFSHFT:',IFSHFT
  write(u6,*) '     IFHDIA:',IFHDIA
  write(u6,*) '     IFHCOM:',IFHCOM
  write(u6,*) '     IFSO  :',IFSO
  write(u6,*) '     NATO  :',NATO
  write(u6,*) '     IFTRD1:',IFTRD1
  write(u6,*) '     IFTRD2:',IFTRD2
  write(u6,*) '     IFTDM :',IFTDM
  write(u6,*) '     RFPERT:',RFPERT
  write(u6,*) '     TOFILE:',ToFile
  write(u6,*) '     PRXVR :',PRXVR
  write(u6,*) '     PRXVE :',PRXVE
  write(u6,*) '     PRXVS :',PRXVS
  write(u6,*) '     PRMER :',PRMER
  write(u6,*) '     PRMEE :',PRMEE
  write(u6,*) '     PRMES :',PRMES
  write(u6,*) '     IFGCAL:',IFGCAL
  write(u6,*) '     IFXCAL:',IFXCAL
  write(u6,*) '     IFMCAL:',IFMCAL
  write(u6,*) '     HOP:',HOP
  write(u6,*) '     TRACK:',TRACK
  write(u6,*) '     ONLY_OVERLAPS:',ONLY_OVERLAPS
  write(u6,*) '     IfDCpl:',IfDCpl
  write(u6,*) '     IFCURD:',IFCURD
  write(u6,*) '     Do_TMOM:',Do_TMOM
  write(u6,*) '     Do_SK:',Do_SK
  write(u6,*) '     L_Eff:',L_Eff
  write(u6,*) '     CD:',DOCD
  write(u6,*) '     Force_NON_AO_TDM:',Force_NON_AO_TDM
end if

! DEFAULT WAVE FUNCTION TYPE:
WFTYPE = 'GENERAL'
if (IPGLOB > 3) write(u6,*) ' ***** INIT ENDS **********'

end subroutine INIT_RASSI
