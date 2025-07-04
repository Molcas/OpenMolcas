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
      SUBROUTINE INIT_RASSI()

      use symmetry_Info, only: symmetry_info_get
      use rassi_aux, only: ipglob
#ifndef _DMRG_
      use rasscf_global, only: doDMRG
#endif
      use cntrl, only: SONTOSTATES, SONATNSTATE, SODIAGNSTATE
      use Cntrl, only: NJOB, NSTATE, NPROP, NSOPR, CITHR,
     &                 NSOThr_Prt, SOThr_Prt, MXPROP, PRSXY, PRDIPVEC,
     &                 PRORB, PRTRA, PRCI, CIH5, IFHAM, IFEJOB, IFSHFT,
     &                 IFHDIA, IFHEXT, IFHEFF, IFHCOM, HAVE_HEFF,
     &                 HAVE_DIAG, NOHAM, IFSO, IFNTO, NATO, BINA,
     &                 IFTRD1, IFTRD2, IFTDM, RFPert, ToFile, PRXVR,
     &                 PRXVS, PRMER, PRMEE, PRMES, IFGCAL, EPRThr,
     &                 EPRAThr, IFXCAL, IFMCAL, HOP, TRACK,
     &                 ONLY_OVERLAPS, DIPR, OSThr_DipR, QIPR,
     &                 OSThr_QIPR, QIALL, DYSO, DYSEXPORT, TDYS, OCAN,
     &                 DCHS, DCHO, DO_TMOM, PRRAW, PRWEIGHT, TOLERANCE,
     &                 REDUCELOOP, LOOPDIVIDE, LOOPMAX, TMGR_Thrs,
     &                 Do_SK, Do_Pol,
     &                 L_Eff, DoCD, RSThr, RSPR, FORCE_NON_AO_TDM,
     &                 IFDCPL, LPRPR, LHAMI, IFATCALSA, IFGTCALSA,
     &                 IFGTSHSA, IFACAL, IFACALFC, IFACALSD, NOSO,
     &                 IFCURD, IFArgU, NrNATO, NBINA, TDIPMIN, JBNAME,
     &                 PNAME, PTYPE, SOPRNM, SOPRTP, PRXVE, MINAME
      use cntrl, only: LuOne, FnOne, LuOrd, FnOrd, LuIph,
     &                 LuExc, FnExc, LuMck, LuTOM, FnTOM,
     &                 LuEig, FnEig
      use rassi_data, only: WFTYPE
      use hfc_logical, only: MAG_X2C


      IMPLICIT None
      Character(LEN=256) STRING
      Logical FoundTwoEls,DoCholesky
      Integer I, IPROP


* Initialise doDMRG if compiled without QCMaquis
#ifndef _DMRG_
      DoDMRG = .false.
#endif

      Call symmetry_info_get()

C UNIT NUMBERS AND NAMES
      LUONE=2
      FNONE='ONEINT'
      LUORD=30
      FNORD='ORDINT'
      LUIPH=15
      LUEXC=22
      FNEXC='ANNI'
      LUMCK=33
      LuToM=26
      FnToM='TOFILE'
      LuEig=27
      FnEig='EIGV'
      JBNAME(:)='UNDEFINE'
      DO  I=1,SIZE(JBNAME)
        WRITE(MINAME(I),'(''MCK'',I3.3)') I
      END DO
      IF(IPGLOB.GT.3) THEN
        WRITE(6,*)' Unit numbers and names:'
        WRITE(6,'(1x,I8,5x,A8)')LUONE,FNONE
        WRITE(6,'(1x,I8,5x,A8)')LUORD,FNORD
        WRITE(6,'(1x,I8,5x,A8)')LUEXC,FNEXC
      END IF

      IF(IPGLOB.GT.3) WRITE(6,*)' OPENING ',FNEXC
      CALL DANAME(LUEXC,FNEXC)


C NR OF JOBIPHS AND STATES:
      NJOB=0
      NSTATE=0
      IF(IPGLOB.GT.3) THEN
       WRITE(6,*)' INITIAL DEFAULT VALUES:'
       WRITE(6,'(1X,A,I4)')'  NJOB:',NJOB
       WRITE(6,'(1X,A,I4)')'NSTATE:',NSTATE
      END IF
C
C NR OF OPERATORS FOR WHICH MATRIX ELEMENTS ARE TO BE CALCULATED:
      NPROP=0

C OPERATORS FOR WHICH MATRIX ELEMENTS OVER SPIN-ORBIT EIGENSTATES
C ARE TO BE COMPUTED.
      NSOPR=0

C DEFAULT THRESHOLD FOR PRINTING CI COEFFICIENTS:
      CITHR=0.05d0

C DEFAULT THRESHOLD FOR PRINTING TRANSITION DIPOLE VECTORS
      TDIPMIN=1.0D-4

C DEFAULT THRESHOLD AND MAX NUMBER OF SO-HAMILTONIAN
C MATRIX ELEMENTS TO PRINT:
      NSOTHR_PRT=0
      SOTHR_PRT=-1.0D0

C SET LABELS TO UNDEFINED
      DO IPROP = 1, MXPROP
         PNAME(IPROP) ='UNDEF.  '
         PTYPE(IPROP) ='UNDEF.  '
         SOPRNM(IPROP)='UNDEF.  '
         SOPRTP(IPROP)='UNDEF.  '
      END DO

C DEFAULT FLAGS:
      PRSXY=.FALSE.
      PRDIPVEC=.FALSE.
      PRORB=.FALSE.
      PRTRA=.FALSE.
      PRCI=.FALSE.
      CIH5=.FALSE.
      IFHAM=.FALSE.
      IFEJOB=.FALSE.
      IFSHFT=.FALSE.
      IFHDIA=.FALSE.
      IFHEXT=.FALSE.
      IFHEFF=.FALSE.
      IFHCOM=.FALSE.
      HAVE_HEFF=.FALSE.
      HAVE_DIAG=.FALSE.
      NOHAM=.FALSE.
      IFSO=.FALSE.
      IFNTO=.FALSE.
      NATO=.FALSE.
      BINA=.FALSE.
      IFTRD1=.FALSE.
      IFTRD2=.FALSE.
      IFTDM=.FALSE.
      RFPERT=.FALSE.
      ToFile=.false.
      PRXVR=.FALSE.
      PRXVE=.FALSE.
      PRXVS=.FALSE.
      PRMER=.FALSE.
      PRMEE=.FALSE.
      PRMES=.FALSE.
      IFGCAL=.FALSE.
      EPRTHR=0.0D0
      EPRATHR=0.0D0
      IFXCAL=.FALSE.
      IFMCAL=.FALSE.
      HOP=.FALSE.
      TRACK=.FALSE.
      ONLY_OVERLAPS=.FALSE.
* Intensities
      DIPR=.FALSE.
      OSTHR_DIPR = 0.0D0
      QIPR=.FALSE.
      OSTHR_QIPR = 0.0D0
      QIALL=.FALSE.
      DYSO=.FALSE.
      DYSEXPORT=.FALSE.
      TDYS=.FALSE.
      OCAN=1
      DCHS=.FALSE.
      DCHO=1
* Exact operator
      Do_TMOM=.FALSE.
      PRRAW=.FALSE.
      PRWEIGHT=.FALSE.
      TOLERANCE=0.1D0
      REDUCELOOP=.FALSE.
      LOOPDIVIDE=0
      LOOPMAX=-1
      TMGr_thrs=-1.0d0
      Do_SK  =.FALSE.
      Do_Pol  =.FALSE.
      L_Eff=5
C CD - velocity and mixed gauge
      DOCD = .FALSE.
      RSTHR = 0.0D0
      RSPR=.FALSE.
C Force that TDMs are not stored in the AO basis.
      Force_NON_AO_TDM=.False.
      CALL GETENVF('MOLCAS_FORCE_NON_AO_TDM',STRING)
      If (STRING.eq.'ON') Force_NON_AO_TDM=.True.
cnf
      IfDCpl = .False.
cnf

C tjd- BMII: Print out spin-orbit properties to files
      LPRPR=.FALSE.
      LHAMI=.FALSE.
c Feng: test control
      MAG_X2C=.FALSE.

C K. Sharkas  BEG
      IFATCALSA=.FALSE.
      IFGTCALSA=.FALSE.
      IFGTSHSA=.FALSE.
C K. Sharkas  END

c BP - Hyperfine tensor and SONATORB initialization
c RF - SO-NTO initialization
      IFACAL=.FALSE.
      IFACALFC=.TRUE.
      IFACALSD=.TRUE.

      NOSO=.FALSE.
      SONATNSTATE=0
      SODIAGNSTATE=0

      SONTOSTATES=0

      IFCURD=.FALSE.
      IFARGU=.FALSE.


* Nr of states for which natural orbitals will be computed:
      NRNATO=0
* Nr of state pairs for computing bi-natural orbitals:
      NBINA=0

* Check if two-electron integrals are available:
      Call f_Inquire('ORDINT',FoundTwoEls)
      Call DecideOnCholesky(DoCholesky)
      If (FoundTwoEls .or. DoCholesky) IFHAM=.True.

      IF(IPGLOB.GE.4) THEN
        WRITE(6,*)'Initial default flags are:'
        WRITE(6,*)'     PRSXY :',PRSXY
        WRITE(6,*)'     PRORB :',PRORB
        WRITE(6,*)'     PRTRA :',PRTRA
        WRITE(6,*)'     PRCI  :',PRCI
        WRITE(6,*)'     IFHAM :',IFHAM
        WRITE(6,*)'     IFHEXT:',IFHEXT
        WRITE(6,*)'     IFHEFF:',IFHEFF
        WRITE(6,*)'     IFEJOB:',IFEJOB
        WRITE(6,*)'     IFSHFT:',IFSHFT
        WRITE(6,*)'     IFHDIA:',IFHDIA
        WRITE(6,*)'     IFHCOM:',IFHCOM
        WRITE(6,*)'     IFSO  :',IFSO
        WRITE(6,*)'     NATO  :',NATO
        WRITE(6,*)'     IFTRD1:',IFTRD1
        WRITE(6,*)'     IFTRD2:',IFTRD2
        WRITE(6,*)'     IFTDM :',IFTDM
        WRITE(6,*)'     RFPERT:',RFPERT
        WRITE(6,*)'     TOFILE:',ToFile
        WRITE(6,*)'     PRXVR :',PRXVR
        WRITE(6,*)'     PRXVE :',PRXVE
        WRITE(6,*)'     PRXVS :',PRXVS
        WRITE(6,*)'     PRMER :',PRMER
        WRITE(6,*)'     PRMEE :',PRMEE
        WRITE(6,*)'     PRMES :',PRMES
        WRITE(6,*)'     IFGCAL:',IFGCAL
        WRITE(6,*)'     IFXCAL:',IFXCAL
        WRITE(6,*)'     IFMCAL:',IFMCAL
        WRITE(6,*)'     HOP:',HOP
        WRITE(6,*)'     TRACK:',TRACK
        WRITE(6,*)'     ONLY_OVERLAPS:',ONLY_OVERLAPS
        WRITE(6,*)'     IfDCpl:',IfDCpl
        WRITE(6,*)'     IFCURD:',IFCURD
        WRITE(6,*)'     Do_TMOM:',Do_TMOM
        WRITE(6,*)'     Do_SK:',Do_SK
        WRITE(6,*)'     L_Eff:',L_Eff
        WRITE(6,*)'     CD:',DOCD
        WRITE(6,*)'     Force_NON_AO_TDM:',Force_NON_AO_TDM
      END IF

C DEFAULT WAVE FUNCTION TYPE:
      WFTYPE='GENERAL '
      IF(IPGLOB.GT.3) WRITE(6,*)' ***** INIT ENDS **********'

      END SUBROUTINE INIT_RASSI
