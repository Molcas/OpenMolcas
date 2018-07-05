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
      Subroutine Proc_InpX(DSCF,Info,lOPTO,iRc)

! module dependencies
#ifdef module_DMRG
!     use molcas_dmrg_interface !stknecht: Maquis-DMRG program
#endif

      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
#include "rasdim.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "gas.fh"
#include "rasscf.fh"
#include "input_ras.fh"
#include "splitcas.fh"
#include "bk_approx.fh"
#include "general.fh"
#include "output_ras.fh"
#include "orthonormalize.fh"
      Parameter (ROUTINE='READIN  ')
#include "casvb.fh"
#include "pamint.fh"
* Lucia-stuff:
#include "ciinfo.fh"
#include "csfbas.fh"
#include "spinfo.fh"
#include "lucia_ini.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      character(32) :: prgm
#endif
*
      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      Character*16  OFE_KSDFT
      COMMON  / OFembed_C / OFE_KSDFT
      COMMON  / OFembed_I / ipFMaux, ip_NDSD, l_NDSD
      COMMON  / OFembed_T / ThrFThaw
      COMMON  / OFembed_R1/ Xsigma
      COMMON  / OFembed_R2/ dFMD
*
      Character*180  Line
!      Character*8 NewJobIphName
      logical lExists, RunFile_Exists, RlxRCheck
* Some strange extra logical variables...
      logical lOPTO
      logical DSCF
      logical RF_On
      logical Langevin_On
      logical PCM_On
      Integer ipTemp1,ipTemp2,ipTemp3
* (SVC) added for treatment of alter and supsym
      Dimension iMAlter(8,2)
      Integer IPRGLB_IN, IPRLOC_IN(7)

      Logical DoCholesky,timings,DensityCheck
      Logical DoLocK,Deco
      Logical Estimate,Update
      Integer ALGO,Nscreen
      Real*8  dmpk,ChFracMem
      Logical DBG, exist

      Common /CHLCAS / DoCholesky,ALGO
      COMMON /CHODENSITY/ DensityCheck
      COMMON /CHOTIME / timings
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen
      COMMON /CHOSCREEN/ Estimate,Update
      COMMON /CHOPAR/ ChFracMem

      Integer IScratch(10)
* Label informing on what type of data is available on an INPORB file.
      Character*8 InfoLbl
* Local NBAS_L, NORB_L .. avoid collision with items in common.
      DIMENSION NBAS_L(8),NORB_L(8)
      DIMENSION NFRO_L(8),NISH_L(8),NRS1_L(8),NRS2_L(8)
      DIMENSION NRS3_L(8),NSSH_L(8),NDEL_L(8)
#ifdef _HDF5_
      character(1), allocatable :: typestring(:)
#endif
* TOC on JOBOLD (or JOBIPH)
      DIMENSION IADR19(15)

      Character*180 Get_LN
      External Get_LN
      Real*8   Get_ExFac
      External Get_ExFac
      Character*100 ProgName, Get_SuperName
      External Get_SuperName
      Character*72 ReadStatus
      Character*72 JobTit(mxTit)
      Character*256 myTitle
      Character*256 RealName
      Logical, External :: Is_First_Iter

      INTEGER :: iDNG
      Logical :: DoDMRG,DNG
      Character*8 emiloop
      Character*8 inGeo
      Logical :: MustCopy

      Intrinsic INDEX,NINT,DBLE,SQRT
C...Dongxia note for GAS:
C   No changing about read in orbital information from INPORB yet.

      Call StatusLine('MCPDFT:','Processing Input')

!      DBG = .TRUE.
      DBG = .FALSE.
      DoFaro = .FALSE.

* NN.14 Block DMRG flag
      DoDMRG = .false.

      doGradPDFT = .false.
      doNOGRad = .false.

* Orbital-free embedding
      Do_OFemb=.false.
      KEonly  =.false.
      OFE_first  =.true.
      ipFMaux = -666666
      ip_NDSD = -696969
      l_NDSD = 0
      ThrFThaw = 0.0d0
      dFMD = 0.0d0
      Xsigma=1.0d4

*    SplitCAS related variables declaration  (GLMJ)
      DoSplitCAS= .false.
      NumSplit = .false.
      EnerSplit = .false.
      PerSplit = .false.
      FOrdSplit  = .false.
*    BK type of approximation (GLMJ)
      DoBKAP = .false.

*    GAS flag, means the INPUT was GAS
      iDoGas = .false.

      Call qEnter('Proc_Inp')

      NAlter=0
      iRc=_RC_ALL_IS_WELL_

* Note: During process of keywords, there is either a call to ChkIfKey
* once any input data for that keyword has been processed, to check that
* there are no unrecognized keywords (misspelled? syntactic error?)
* following; or if the keyword should not be followed by any data,
* then for the same purpose, a call of the kind
*       Call SetPos(LUInput,'ATOM',Line,iRc)
*       Call ChkIfKey()
* This is really the only safe way to make warnings for such events
* the way the input is parsed and processed for the moment.
* Start by checking for extraneous input following any of these
* keywords (misspelled or nonsyntactic usage)
      If (KeyCORE) Then
        Call SetPos_m(LUInput,'CORE',Line,iRc)
        Call ChkIfKey_m()
      End If
      If (KeyINPO) Then
        Call SetPos_m(LUInput,'INPO',Line,iRc)
        Call ChkIfKey_m()
      End If
      If (KeyJOBI) Then
        Call SetPos_m(LUInput,'JOBI',Line,iRc)
        Call ChkIfKey_m()
      End If
      If (KeyLUMO) Then
        Call SetPos_m(LUInput,'LUMO',Line,iRc)
        Call ChkIfKey_m()
      End If
      If (KeyCIRE) Then
        Call SetPos_m(LUInput,'CIRE',Line,iRc)
        Call ChkIfKey_m()
      End If
      !ANDREWANDREWANDREW
      KeyCIRE=.TRUE.

* They had to be tested here, since the key flags may be modified later.
*
* ====================================================================
* The outdated INPORB keyword now means just the same as LUMORB
* Will be deleted soon.
      KeyLUMO=KEYLUMO .or. KeyINPO
* ====================================================================
*
* How was the program called?
*PAM 2009  For some particular types of calculations, the input is to be
* disregarded or overridden, as follows (Roland, esp. numerical diff):
      MustCopy=.False.
      If(KeyEXPE) Then
       Call SetPos_m(LUInput,'EXPE',Line,iRc)
       Call ChkIfKey_m()
      Else
        ProgName=Get_SuperName()
        IfVB=0
        If (ProgName(1:6).eq.'mcpdft') Then
* For geometry optimizations use the old CI coefficients.
         If (.Not.Is_First_Iter()) Then
           KeyCIRE=.true.
           KeyFILE=.false.
           MustCopy=.true.
         End If
        Else If (ProgName(1:5).eq.'casvb') Then
         IfVB=2
        Else If (ProgName(1:6).eq.'loprop') Then
         KeyCIRE=.true.
         KeyFILE=.false.
        Else If (ProgName(1:11).eq.'last_energy') Then
         KeyCIRE=.true.
         KeyFILE=.false.
         MustCopy=.true.
        Else If (ProgName(1:18).eq.'numerical_gradient') Then
         KeyCIRE=.true.
         KeyFILE=.false.
         MustCopy=.true.
        End If
      End If
*PAM2009 Also, CIRESTART would normally also imply that orbitals are
*to be taken from JOBIPH or JOBOLD:
      If(.not.KeyEXPE) Then
        If(KeyCIRE) Then
         KeyLUMO=.false.
         KeyINPO=.false.
         KeyJOBI=.true.
        End If
      End If

* Check PRINT command
      If (KeyPRIN) Then
       Call SetPos_m(LUInput,'PRIN',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after PRINT keyword.'
       Read(LUInput,*,End=9910,Err=9920)
     &          LF,IPRGLB_IN,(IPRLOC_IN(I),I=1,7)
       ReadStatus=' O.K. reading info of PRINT keyword.'
       Call SetPrLev_m(LF,IPRGLB_IN,IPRLOC_IN)
       If (IPRLOC(1).GE.DEBUG) Then
        Write(6,*)' PROC_INP: Print levels have been set:'
        Write(6,*)'    Global print level IPRGLB=',IPRGLB
        Write(6,*)'    Local print levels by section:'
        Write(6,*)'             Input section, IPRLOC(1)=',IPRLOC(1)
        Write(6,*)'    Transformation section, IPRLOC(2)=',IPRLOC(2)
        Write(6,*)'                CI section, IPRLOC(3)=',IPRLOC(3)
        Write(6,*)'          Super-CI section, IPRLOC(4)=',IPRLOC(4)
        Write(6,*)'            Output section, IPRLOC(6)=',IPRLOC(6)
        Write(6,*)'          Property section, IPRLOC(7)=',IPRLOC(7)
       End If
       Call ChkIfKey_m()
      End If

* Local print level in this routine:
!AMS
      IPRLEV=IPRLOC(1)
!      IPRLEV=INSANE
* Short, for use with IF-statements at debugging print level:
      DBG=DBG .or. (IPRLEV.GE.DEBUG)

      If (DBG) Then
        Write(6,*)' Trace of input processing in PROC_INP:'
      End If

* ====================================================================
* Check IfVB flag
* If (IFVB.eq.2) -- Then bypass a lot of processing:
      If (DBG) Write(6,*)' Valence-Bond flag IfVB=',IfVB
      if(ifvb.eq.2)then
        if(JOBIPH.le.0) Then
          JOBIPH=IsFreeUnit(15)
          call daname(JOBIPH,"JOBIPH")
        end if
* Special input routine:
        If (DBG) Write(6,*)' IfVB=2, so Call Readin_VB'
        Call Readin_vb_m()
        If (DBG) Write(6,*)' Back from Readin_VB'
        If (DBG) Write(6,*)' Bypass usual input processing, GoTo 100!'
      endif
      if(ifvb.eq.2) GoTo 100

* ====================================================================
      If (DBG) Write(6,*)' Check if VB keyword was given.'
      If (KeyVB) Then
       If (DBG) Write(6,*)' Yes it was!'
*---  Process vb   command --------------------------------------------*
       Call SetPos_m(LUInput,'VB  ',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       If (DBG) Write(6,*)' so Call CVBInp_CVB.'
       call cvbinp_rvb_m(1,LUInput)
       If (DoCholesky) then
        Call WarningMessage(2,'CASVB cannot do Cholesky or RI/DF.')
        Call Quit(_RC_INPUT_ERROR_)
       End If
       If (DBG) Write(6,*)' Set IfVB=1.'
       ifvb=1
       If (DBG) Write(6,*)' Asked for CASVB calc.'
      End If

* ====================================================================
      If (KeyTITL) Then
*---  process TITLE    command ----------------------------------------*
       Call SetPos_m(LUInput,'TITL',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Title(1)=Trim(Get_Ln(LUInput))
       If (DBG) Then
         Write(6,*)' PROC_INP: Title line:'
         Write(6,*) Title(1)
       End If
       Call ChkIfKey_m()
      End If

* ====================================================================
*---  Check if there is any runfile:
      Call F_Inquire('RUNFILE',RunFile_Exists)
      If (DBG) Write(6,*)' Inquire about RUNFILE.'
      IF (RunFile_Exists) Then
       If (DBG) Write(6,*)' Yes, there is one.'
       NSYM=0
       Call qpg_iScalar('nSym',lExists)
       IF (lExists) Then
        Call Get_iScalar('nSym',nSym)
        Call Get_iArray('nBas',nBas,nSym)
        If (DBG) Then
          Write(6,*)' The following information exists on runfile:'
          Write(6,*)' Nr of symmetries, NSYM:',NSYM
          Write(6,*)' Nr of basis functions/symmetry:'
          Write(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
          Call XFlush(6)
        End If
       ELSE
        Call WarningMessage(2,'No symmetry info on runfile.')
        Write(6,*)' There seems to be no information about symmetry'
        Write(6,*)' on the runfile! This is an unexpected error.'
        Call Quit(_RC_IO_ERROR_READ_)
       END IF
      ELSE
       Call WarningMessage(2,'Cannot find runfile.')
       Write(6,*)' PROC_INP: Cannot find RUNFILE. This is an'//
     &           ' unexpected error.'
        Call Quit(_RC_IO_ERROR_READ_)
      END IF

* PAM03: Local variable IORBDATA=0 as long as no explicit orbital specifications have
* been set by the user input.
* iOrbData=0: no orbital space data is specified
*          1: user input has defined orbital spaces
*          2: specifications from runfile, or JOBOLD or JOBIPH file
*          3: specifications from orbital file
      iOrbData=0
      INVEC=0
      iHAVECI=0
* INVEC=0, no source for orbitals (yet)
*       1, CORE command: compute orbitals from scratch.
*       2, read from starting orbitals file in INPORB format.
*       3, take from JOBOLD, or JOBIPH file
* RASSCF will read CMO arrays by READVC call after return from PROC_INP
* Inside readvc, if INVEC is still 0, it will be changed to 4,5, or 6:
*       4, take from an HDF5 file
*       5, take rasscf orbitals from runfile
*       6, take scf orbitals from runfile
*       7, take guessorb orbitals from runfile
* iOverWr: If orbital coefficient arrays (CMO:s) are read from INPORB
* file, then iOverWr=0 means that the arrays are reordered by typeindex.
* iOverWr=1 means do not reorder by typeindex.
      iOverWr=0

*--- Where should orbital be read (if at all)?
      If (DBG) Write(6,*)' Where to read orbitals? '
      StartOrbFile='INPORB'
      If (KeyFILE) Then
       If (DBG) Then
         Write(6,*)' Reading file name for start orbitals.'
       End If
       Call SetPos_m(LUInput,'FILE',Line,iRc)
       Line=Get_Ln(LUInput)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Call ChkIfKey_m()
       If (DBG) Then
         Write(6,*) ' Calling fileorb with filename='
         Write(6,*) Line
       End If
       call fileorb(Line,StartOrbFile)
#ifdef _HDF5_
       if (mh5_is_hdf5(StartOrbFile)) then
         KeyLUMO=.false.
         KeyTYPE=.false.
         KeyCIRE=.false.
         KeyH5OR=.true.
       else
         KeyLUMO=.true.
       end if
#else
       KeyLUMO=.true.
#endif
      End If
      If (DBG) Then
         Write(6,*)' StartOrbFile='
         Write(6,*) StartOrbFile
         call xflush(6)
      End If

* ============================================================================
* The JOBIPH file, following decisions from the Torre Normanna labour camp:
* Default, the file name is 'JOBIPH'.
* However, if keyword IPHNAME was used, then the name was given in input.
* Also, if instead the keyword NEWIPH was given, then a new name will be
* chosen as the first not-already-used name in the sequence
* 'JOBIPH', 'JOBIPH01', 'JOBIPH02', etc.

!Check if JOBOLD exists
      Call F_Inquire('JOBIPH',lExists)
      If(lExists) Then
        Call PrgmTranslate('JOBIPH',RealName,not_sure)
        If (DBG) Then
          Write(6,*)' A JOBIPH file has been found:'
          write(6,*) RealName
        End If
!        JOBOLD=IsFreeUnit(19)
!        Call DANAME(JOBOLD,RealName)
      Else
          Write(LF,*)
          Write(LF,*)'******************************************'
          Write(LF,*)' JOBIPH does not seem to exist,           '
          Write(LF,*)' so the calculation cannot continue.      '
          Write(LF,*)'******************************************'
          Call Abend()
      End If
!      If (DBG) Then
!        Write(6,*)' A fresh JOBPDFT file may be produced.'
!      End If
      IPHNAME='ToBeFoun'
      If (KeyIPHN) Then
        If (DBG) Then
          Write(6,*)' Reading file name for JOBIPH file.'
        End If
        Call SetPos_m(LUInput,'IPHN',Line,iRc)
        If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
        ReadStatus=' Failure reading IPHNAME string.'
        Read(LUInput,*,End=9910,Err=9920) IPHNAME
        ReadStatus=' O.K. after reading IPHNAME string.'
        Call UpCase(IPHNAME)
      End If
      IF(IPHNAME.EQ.'ToBeFoun') THEN
!        IPHNAME='JOBPDFT'
      IPHNAME='JOBIPH'
      END IF
      if(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
      end if
      JOBIPH=IsFreeUnit(15)
!      CALL DANAME(JOBIPH,IPHNAME)
      CALL DANAME(JOBIPH,"JOBIPH")

!      Call F_Inquire('JOBOLD',lExists)
!      If(lExists.and..not.MustCopy) Then
!        If (DBG) Then
!          Write(6,*)' A JOBOLD file has been given.'
!        End If
!      Else !JOBOLD does not exist
!        If(MustCopy) Then
!          If (DBG) Then
!            Write(6,*)' Must overwrite JOBOLD with JOBIPH for grad calc'
!          End If
!        Else
!          If (DBG) Then
!            Write(6,*)' No JOBOLD file, checking for JOBIPH.'
!          End If
!        End If
!        Call F_Inquire('JOBIPH',lExists2)
!        If(lExists2) then
!          If (DBG) Then
!            Write(6,*)' JOBIPH found, copy to JOBOLD.'
!          End If
!          Call PrgmTranslate('JOBIPH',RealName,not_sure)
!          write(*,*) 'realname',RealName
!!          Call fcopy(RealName,'JOBOLD',iErr)
!        Else
!          Write(LF,*)
!          Write(LF,*)'******************************************'
!          Write(LF,*)' Neither JOBIPH or JOBOLD seem to exist,  '
!          Write(LF,*)' so the calculation cannot continue.      '
!          Write(LF,*)'******************************************'
!          Call Abend()
!        End If
!      End If
!      If (DBG) Then
!        Write(6,*)' A fresh JOBIPH file will be used.'
!      End If
!      IPHNAME='ToBeFoun'
!      If (KeyIPHN) Then
!        If (DBG) Then
!          Write(6,*)' Reading file name for JOBIPH file.'
!        End If
!        Call SetPos_m(LUInput,'IPHN',Line,iRc)
!        If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
!        ReadStatus=' Failure reading IPHNAME string.'
!        Read(LUInput,*,End=9910,Err=9920) IPHNAME
!        ReadStatus=' O.K. after reading IPHNAME string.'
!        Call UpCase(IPHNAME)
!      End If
!* The JOBIPH file:
!      IF(IPHNAME.EQ.'ToBeFoun') THEN
!* Choose a new jobiph name.
!        IPHNAME='JOBPDFT'
!      END IF
!      if(JOBIPH.gt.0) Then
!        Call DaClos(JOBIPH)
!        JOBIPH=-1
!      end if
!      JOBIPH=IsFreeUnit(15)
!      CALL DANAME(JOBIPH,IPHNAME)



!      If (DBG) Write(6,*)' Present name of JOBIPH file is ',IPHNAME
!      Call F_Inquire('JOBIPH',lExists)
!      IF (DBG.and.lExists) Write(6,*) 'JOBIPH exists; may be read.'
*---  Process NEWIPH command (P A Malmqvist Sep 06)------------------------*
!      KeyNEWI=.TRUE.
!      If (KeyNEWI) Then
!       If (DBG) Then
!        Write(6,*)' A fresh JOBIPH file will be used.'
!       End If
!       IPHNAME='ToBeFoun'
!      End If
*---  process IPHNAME command (P A Malmqvist Sep 06)------------------------*
!      If (KeyIPHN) Then
!       If (DBG) Then
!        Write(6,*)' Reading file name for JOBIPH file.'
!       End If
!       Call SetPos_m(LUInput,'IPHN',Line,iRc)
!       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
!       ReadStatus=' Failure reading IPHNAME string.'
!       Read(LUInput,*,End=9910,Err=9920) IPHNAME
!       ReadStatus=' O.K. after reading IPHNAME string.'
!       Call UpCase(IPHNAME)
!      End If
* The JOBIPH file:
!      IF(IPHNAME.EQ.'ToBeFoun') THEN
* Choose a new jobiph name.
!       Call F_Inquire('JOBIPH',lExists)
!       IF (.not.lExists) THEN
!        IPHNAME='JOBIPH'
!       ELSE
!        DO I=1,99
!         write(NewJobIphName,'(a6,i2.2)') 'JOBIPH',I
!         Call F_Inquire(NewJobIphName,lExists)
* Choose first non-existent file
!         IF(.not.lExists) goto 12
!        END DO
!        Write(LF,*)
!        Write(LF,*)'******************************************'
!        Write(LF,*)' Sorry, all possible JOBIPH names in this '
!        Write(LF,*)' directory already in use. Program stops. '
!        Write(LF,*)'******************************************'
!        Call Abend()
!  12    Continue
!        IPHNAME=NewJobIphName
!       END IF
!      END IF
      !If (DBG) Then
!        Write(6,*) ' Name of JOBIPH file is'
!        Write(6,*) IPHNAME
      !End If
* Finally, we know the name. Open jobiph file. If another file is already
* opened with this identifier close it.
!      if(JOBIPH.gt.0) Then
!        Call DaClos(JOBIPH)
!        JOBIPH=-1
!      end if
!      JOBIPH=IsFreeUnit(15)
!      CALL DANAME(JOBIPH,IPHNAME)

* ========================================================================
*  If orbital files should be produced, several keywords are relevant:
*---  Process OUTPRINT command -----
      If (KeyOUTP) Then
       If (DBG) Write(6,*) ' OUTPRINT command was given:'
       Call SetPos_m(LUInput,'OUTP',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading OUTPRINT line.'
       Read(LUInput,*,End=9910,Err=9920) Line
       ReadStatus=' O.K. after reading OUTPRINT line.'
       Call UpCase(Line)
       If ( Index(Line,'VERY').ne.0 ) OutFmt1='NOTHING '
       If ( Index(Line,'BRIE').ne.0 ) OutFmt1='FEW     '
       If ( Index(Line,'LONG').ne.0 ) OutFmt1='ALL     '
       If ( Index(Line,'NOTH').ne.0 ) OutFmt1='NOTHING '
       If ( Index(Line,'FEW' ).ne.0 ) OutFmt1='FEW     '
       If ( Index(Line,'NOCO').ne.0 ) OutFmt1='NOCORE  '
       If ( Index(Line,'ALL' ).ne.0 ) OutFmt1='ALL     '
       If ( Index(Line,'COMP').ne.0 ) OutFmt2='COMPACT '
       If ( Index(Line,'FULL').ne.0 ) OutFmt2='FULL    '
       If (OutFmt1.eq.'DEFAULT '.and.OutFmt2.eq.'DEFAULT ') Then
        Call WarningMessage(1,'Error in ''OUTP'' command?')
        Write(LF,*)' Input line is:'
        Write(LF,*) Line
        Write(LF,*)' Did not understand ''OUTP'' input. Ignored.'
       End If
      End If
*---  Process PROR (Print levels for orbitals) command --- (new!) -----*
      If (KeyPROR) Then
       If (DBG) Write(6,*) ' PRORB command was given:'
       Call SetPos_m(LUInput,'PROR',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading PRORB input.'
       Read(LUInput,*,End=9910,Err=9920) PRETHR,PRO
       ReadStatus=' O.K. after reading PRORB input.'
       PROTHR=MAX(0.0D0,PRO)
      End If
*---  Process ORBL keyword: Orbital listing
      If (KeyORBL) Then
       If (DBG) Write(6,*) ' ORBL (Orbital listing):'
       Call SetPos_m(LUInput,'ORBL',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading ORBL input.'
       Read(LUInput,*,End=9910,Err=9920) Line
       ReadStatus=' O.K. after reading ORBL input.'
       Call UpCase(Line)
       If ( Index(Line,'VERY').ne.0 ) OutFmt1='NOTHING '
       If ( Index(Line,'BRIE').ne.0 ) OutFmt1='FEW     '
       If ( Index(Line,'LONG').ne.0 ) OutFmt1='ALL     '
       If ( Index(Line,'NOTH').ne.0 ) OutFmt1='NOTHING '
       If ( Index(Line,'FEW' ).ne.0 ) OutFmt1='FEW     '
       If ( Index(Line,'NOCO').ne.0 ) OutFmt1='NOCORE  '
       If ( Index(Line,'ALL' ).ne.0 ) OutFmt1='ALL     '
       If (OutFmt1.eq.'DEFAULT ') Then
        Call WarningMessage(1,'Error in ''ORBL'' command?')
        Write(LF,*)' Input line is:'
        Write(LF,*) Line
        Write(LF,*)' Did not understand ''OUTL'' input. Ignored.'
       End If
      End If
*---  Process ORBA keyword: Orbital Appearance
      If (KeyORBA) Then
       If (DBG) Write(6,*) ' ORBA (Orbital appearance):'
       Call SetPos_m(LUInput,'ORBA',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading ORBA input.'
       Read(LUInput,*,End=9910,Err=9920) Line
       ReadStatus=' O.K. after reading ORBA input.'
       Call UpCase(Line)
       If ( Index(Line,'COMP').ne.0 ) OutFmt2='COMPACT '
       If ( Index(Line,'FULL').ne.0 ) OutFmt2='FULL    '
       If (OutFmt2.eq.'DEFAULT ') Then
        Call WarningMessage(1,'Error in ''OUTA'' command?')
        Write(LF,*)' Input line is:'
        Write(LF,*) Line
        Write(LF,*)' Did not understand ''OUTA'' input. Ignored.'
       End If
      End If
*---  Process MAXO keyword: Max nr of state-specific orbital files produced
      If (KeyMAXO) Then
       If (DBG) Write(6,*) ' MAXORB command:'
       Call SetPos_m(LUInput,'MAXO',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading MAXO input.'
       Read(LUInput,*,End=9910,Err=9920) MAXORBOUT
       ReadStatus=' O.K. after reading MAXO input.'
      End If
      If (DBG) Then
        Write(6,*) ' Orbital print levels are'
        Write(6,*) '   for energy, PRETHR=',PRETHR
        Write(6,*) '   for occup , PROTHR=',PROTHR
        Write(6,*) ' Orbital listing flag is: ',OutFmt1
        Write(6,*) ' Orbital appearance flag: ',OutFmt2
         Write(6,*) ' Max nr of state-specific orbital files is ',
     &                MAXORBOUT
      End If
* ========================================================================
*---  Process OUTORBITALS command (Which kind of orbitals)---------*
      If (KeyOUTO) Then
       Call SetPos_m(LUInput,'OUTO',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading OUTO input.'
       Read(LUInput,*,End=9910,Err=9920) Line
       ReadStatus=' O.K. after reading OUTO input.'
       Call UpCase(Line)
       iOrbTyp=0
       If ( Index(Line,'AVER').ne.0 ) iOrbTyp=1
       If ( Index(Line,'CANO').ne.0 ) iOrbTyp=2
       If ( Index(Line,'NATU').ne.0 ) iOrbTyp=3
       If ( Index(Line,'SPIN').ne.0 ) iOrbTyp=4
       If(iOrbTyp.eq.0) then
         Write(LF,*)' The line after keyword ''OUTORBITALS'' is'
         Write(LF,*)' not understood. That line begins:'
         Write(LF,'(1x,a60)') line(1:60)
         Write(LF,*)' This input is IGNORED.'
       End If
       If (iOrbTyp.eq.2) IPT2=1
       If (iOrbTyp.eq.3 .or. iOrbTyp.eq.4) Then
         ReadStatus=' Failure reading nOrbRoot after OUTO keyword.'
         Read(LUInput,*,End=9910,Err=9920) nOrbRoot
         ReadStatus=' O.K. after reading nOrbRoot after OUTO keyword.'
       End If
      End If
       If (DBG) Then
         Write(6,*) ' OUTORBITALS command specified orbital type ',
     &                iOrbTyp
         IF(iOrbTyp.eq.1) Write(6,*) ' Meaning: Average'
         IF(iOrbTyp.eq.2) Write(6,*) ' Meaning: Canonical'
         IF(iOrbTyp.eq.3) Write(6,*) ' Meaning: Natural'
         IF(iOrbTyp.eq.4) Write(6,*) ' Meaning: Spin orbitals.'
         If (iOrbTyp.eq.3 .or. iOrbTyp.eq.4) Then
          Write(6,*) ' Max state for printing this orbital type ',
     &                nOrbRoot
         End If
       End If
*---  Process ORDER command (SVC Feb 06)----------------------------------*
      If (KeyORDE) Then
       If (DBG) Write(6,*) ' ORDER command was used.'
       Call SetPos_m(LUInput,'ORDE',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure after reading ORDER keyword.'
       Read(LUInput,*,End=9910,Err=9920) IFORDE
       ReadStatus=' O.K. after reading ORDER keyword.'
       IORDEM=1
       If (DBG) Then
        Write(6,*) ' IFORDE, IORDEM=',IFORDE,IORDEM
       End If
      End If
*---  Process PRSP command --------------------------------------------*
      If (KeyPRSP) Then
       If (DBG) Write(6,*) ' PRSPIN command was used.'
       ISPDEN=1
      End If
* ========================================================================
*  If ORBONLY keyword was used, then the JOBIPH file should be used
* only to produce orbital files, then the program stops.
*---  Process ORBO command -( new! generate orbitals from jobiph, only)*
!      If (KeyORBO) Then
!       If (DBG) Write(6,*) ' ORBOnly command was used.'
!       iOrbOnly=1
!       Call OrbFiles_m(JOBIPH,IPRLEV)
!       If(JOBIPH.gt.0) Then
!         Call DaClos(JOBIPH)
!         JOBIPH=-1
!       End If
!* Nothing more to be done, so return.
!       iReturn=_RC_ALL_IS_WELL_
!       Call xQuit(iReturn)
!      End If
*
*---  Process ALTEr command (G. Ghigo Sep 03)--------------------------*
!      If (KeyALTE) Then
!       If (DBG) Write(6,*) ' ALTER command has been used.'
!       Call SetPos_m(LUInput,'ALTE',Line,iRc)
!       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
!       ReadStatus=' Failure after reading ALTER keyword.'
!       Read(LUInput,*,End=9910,Err=9920) NAlter
!       ReadStatus=' O.K. after reading ALTER keyword.'
!       If ( NAlter.gt.8 ) Then
!        Write(LF,*)
!        Call WarningMessage(2,'Alter: too many orbital pairs.')
!        Write(LF,*) ' ************* ERROR **************'
!        Write(LF,*) ' ALTEr: Too many pairs of orbitals '
!        Write(LF,*) ' to exchange (max 8).              '
!        Write(LF,*) ' **********************************'
!        Call Abend()
!       End If
!       Do iAlter=1,NAlter
!        ReadStatus=' Failure reading data after ALTER keyword.'
!        Read(LUInput,*,End=9910,Err=9920) (MAlter(iAlter,i),i=1,3)
!        ReadStatus=' O.K. after reading data after ALTER keyword.'
!       END DO
!* (SVC) get absolute orbital values for the alterations so that
!* iMAlter is symmetry independent
!       If (DBG) Write(6,*)' ''Absolute'' iMAlter indices:'
!       Do iAlter=1,NAlter
!        iEnd=0
!        iStart=1
!        Do iSym=1,MAlter(iAlter,1)
!          iStart=iEnd+1
!          iEnd=iEnd+nBas(iSym)
!        End Do
!        iMAlter(iAlter,1)=MAlter(iAlter,2)+iStart-1
!        iMAlter(iAlter,2)=MAlter(iAlter,3)+iStart-1
!        If (DBG) Write(6,'(1x,2I5)')
!     &                          iMAlter(iAlter,1),iMAlter(iAlter,2)
!       End Do
!      End If
*---  Process ATOM command (P A Malmqvist Apr 07)----------------------*
!      If(KeyATOM) Then
!       PURIFY='ATOM    '
!       ISUPSM=1
!       Call SetPos_m(LUInput,'ATOM',Line,iRc)
!       Call ChkIfKey_m()
!      End If
!*---  Process LINEAR command (P A Malmqvist Apr 05)----------------------*
!      If(KeyLINE) Then
!       PURIFY='LINEAR'
!       ISUPSM=1
!       Call SetPos_m(LUInput,'LINE',Line,iRc)
!       Call ChkIfKey_m()
!      End If
      If (DBG) Write(6,*) ' Purify=',PURIFY

*---  process KSDF command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if KSDFT was requested.'
      If (KeyKSDF) Then
       If (DBG) Write(6,*) ' KSDFT command was given.'
       PamGen=.False.
       PamGen1=.False.
!AMS       DFTFOCK='CAS '
       DFTFOCK='ROKS'
       Call SetPos_m(LUInput,'KSDF',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
!       Read(LUInput,*,End=9910,Err=9920) Line
!       Call UpCase(Line)
!       If (Line(1:4).eq.'ROKS') DFTFOCK='ROKS'
!       If (Line(1:6).eq.'CASDFT') DFTFOCK='DIFF'
       Read(LUInput,*,End=9910,Err=9920) Line
       KSDFT=Line(1:16)
       Call UpCase(KSDFT)
       If (IPRLOC(1).GE.DEBUG) Then
         If(KSDFT(1:5).eq.'TLSDA')  !GLM
     &     write(6,*) ' TLSDA functional aka LSDA for MCPDFT'
         If(KSDFT(1:5).eq.'TBLYP')
     &     write(6,*) ' TBLYP functional aka BLYP for MCPDFT'
         If(KSDFT(1:4).eq.'TPBE')
     &     write(6,*) ' TPBE functional aka PBE for MCPDFT'
         If(KSDFT(1:5).eq.'FTPBE')
     &     write(6,*) ' FTPBE functional aka PBE for MCPDFT'
         If(KSDFT(1:6).eq.'FTBLYP')
     &     write(6,*) ' FTBLYP functional aka BLYP for MCPDFT'
         If(KSDFT(1:6).eq.'FTLSDA')
     &     write(6,*) ' FTLSDA functional aka LSDA for MCPDFT'
       End if
CGG Calibration of A, B, C, and D coefficients in SG's NewFunctional 1
       If ( KSDFT(1:4).eq.'NEWF') Then
         ReadStatus=' Failure reading data following KSDF=NEWF.'
         Read(LUInput,*,End=9910,Err=9920)
     &                                       Acoef,Bcoef,Ccoef,Dcoef
         ReadStatus=' O.K. after reading data following KSDF=NEWF.'
       End If
CGG This part will be removed. (PAM 2009: What on earth does he mean??)
      ExFac=Get_ExFac(KSDFT)
*******
*
* Read numbers, and coefficients for rasscf potential calculations:
* nPAM  - number of potentials
* ipPAM - list of potentials
* CPAM  - coeffcients of potentials
* PamGen - switch to generate grid of Rho, grad ....., ......
*
*******
       If ( KSDFT(1:3).eq.'PAM') Then
        If ( KSDFT(4:4).eq.'G') PamGen =.True.
        If ( KSDFT(4:4).eq.'G') PamGen1=.False.
        call dcopy_(nPAMintg,0.0d0,0,CPAM,1)
        ReadStatus=' Failure reading data following KSDF=PAM.'
        Read(LUInput,*,End=9910,Err=9920) nPAM
        ReadStatus=' O.K. after reading data following KSDF=PAM.'
*        Write(LF,*) ' Number included exponent in PAM=',nPAM
        Do iPAM=1,nPAM
          ReadStatus=' Failure reading data following KSDF=PAM.'
          Read(LUInput,*,End=9910,Err=9920) Line
          ReadStatus=' O.K.after reading data following KSDF=PAM.'
          Call RdPAM_m(Line,ipPAM(iPAM),CPAM(iPAM))
        End Do
       End If
       Call ChkIfKey_m()
       Else
        Call WarningMessage(2,'No KSDFT functional specified')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' KSDFT functional type must be     '
        Write(LF,*) ' specified for MCPDFT calculations '
        Write(LF,*) ' **********************************'
        Call Abend()
      End If
*---  Process CION command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if CIONLY case.'
      If (KeyCION) Then
       If (DBG) Write(6,*) ' CIONLY keyword was used.'
       iCIonly=1
       Call SetPos_m(LUInput,'CION',Line,iRc)
       Call ChkIfKey_m()
      End If
*---  Process RFPE command ----- (new!) -------------------------------*
      If(KeyRFPE) Then
       If (DBG) Then
        Write(6,*) ' RFPERT (Response Field Perturbation)'
       End If
       RFpert=.true.
       Call SetPos_m(LUInput,'RFPE',Line,iRc)
       Call ChkIfKey_m()
      End If
*---  Process NONE command --( non-equilibrium reaction field )--------*
      If(KeyNONE) Then
       If (DBG) Write(6,*) ' Non-equilibrium response'
       NonEq = .True.
       Call SetPos_m(LUInput,'NONE',Line,iRc)
       Call ChkIfKey_m()
      End If
*---  Process RFRO command --------------------------------------------*
      If(KeyRFRO) Then
       Call SetPos_m(LUInput,'RFRO',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading IPCMROOT after RFROOT keyword.'
       Read(LUInput,*,End=9910,Err=9920) IPCMROOT
       ReadStatus=' O.K. reading IPCMROOT after RFROOT keyword.'
*
*      Check that the root value is not changed explicitly by input
*
       jPCMRoot=iPCMRoot
       Call Qpg_iScalar('RF0CASSCF root',Exist)
       If (Exist) Then
          CALL Get_iScalar('RF0CASSCF root',jPCMRoot)
       Else
          CALL Put_iScalar('RF0CASSCF root',iPCMRoot)
       End If
*
       If (jPCMRoot.ne.iPCMRoot) Then
*         Write (*,*) 'iPCMRoot changed by explicitly by input.'
*
*         Value changed explicitly by input. Accept the new value.
*
          CALL Put_iScalar('RF CASSCF root',iPCMRoot)
          CALL Put_iScalar('RF0CASSCF root',iPCMRoot)
          Call Qpg_dArray("RF CASSCF Vector",Exist,mConf)
          If (Exist) Then
             Call Allocate_Work(ipRF,mConf)
             Call FZero(Work(ipRF),mConf)
             Call Put_dArray("RF CASSCF Vector",Work(ipRF),mConf)
             Call Free_Work(ipRF)
          End If
*
       Else
*
*         Value not explicitly changed by input. See if it is
*         changed by the RunFile, if it exists there.
*
          Call Qpg_iScalar('RF CASSCF root',Exist)
          If (Exist) Then
             CALL Get_iScalar('RF CASSCF root',iPCMRoot)
          Else
             CALL Put_iScalar('RF CASSCF root',iPCMRoot)
          End If
*
       End If
*
       If (DBG) Then
        Write(6,*) ' RFROOT command was given.'
        Write(6,*) ' Response field for root number ',IPCMROOT
       End If
       Call ChkIfKey_m()
      End If
*---  Process CIRF command --------------------------------------------*
      If(KeyCIRF) Then
       Call SetPos_m(LUInput,'CIRF',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' O.K. reading after keyword CIRF.'
       Read(LUInput,*,End=9910,Err=9920) ICIRFROOT
       ReadStatus=' O.K. reading after keyword CIRF.'
       Call ChkIfKey_m()
       If (DBG) Then
        Write(6,*) ' CIRFROOT command was given.'
        Write(6,*) ' Response field will follow CISE root: ',ICIRFROOT
       End If
      End If
*---  Process CIRO command --------------------------------------------*
      If (DBG) Write(6,*) ' Check for CIROOTS command.'
      IF(KeyCIRO) Then
       If (DBG) Write(6,*) ' CIROOTS command was given.'
       Call SetPos_m(LUInput,'CIRO',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
        Line=Get_Ln(LUInput)
CBOR.. Modification 001011
        Line(80:80)='0'
        ReadStatus=' Failure reading after CIROOTS keyword.'
        Read(Line,*,Err=9920,End=9920) NROOTS,LROOTS,iall
        ReadStatus=' O.K reading after CIROOTS keyword.'
        If (NROOTS.GT.MXROOT) Then
          WRITE(6,*) "Error: number of roots exceeds maximum"
          WRITE(6,*) "NROOTS = ", NROOTS
          WRITE(6,*) "MXROOT = ", MXROOT
          CALL AbEnd()
        End If
        if(iall.eq.1) then
         Do i=1,NROOTS
          iroot(i)=i
          WEIGHT(i)=1.d0/DBLE(NROOTS)
         END DO
        Else
CBOR.. End modification 001011
         Line=Get_Ln(LUInput)
         ReadStatus=' Failure reading after CIROOTS keyword.'
         Read(Line,*,Err=9920,End=9920) (IROOT(I),I=1,NROOTS)
         ReadStatus=' O.K.after CIROOTS keyword.'
         Call dCopy_(mxRoot,0.0D0,0,WEIGHT,1)
         If ( NROOTS.eq.1 ) then
           WEIGHT(1)=1.0D0
         Else
           Call GetMem('Temp1','Allo','Inte',ipTemp1,NROOTS)
           Line=Get_Ln(LUInput)
           ReadStatus=' Failure reading after CIROOTS keyword.'
           Read(Line,*,Err=9920) (iWork(ipTemp1+i-1),i=1,NROOTS)
           ReadStatus=' O.K.after CIROOTS keyword.'
           iSum=0
           Do i=1,NROOTS
              iSum=iSum+iWork(ipTemp1+i-1)
           End Do
           Do i=1,NROOTS
              WEIGHT(i)=DBLE(iWork(ipTemp1+i-1))/DBLE(iSum)
           End Do
           Call GetMem('Temp1','Free','Inte',ipTemp1,NROOTS)
         End If
        End If
        If (DBG) Then
         Write(6,*) ' Nr of roots in CI: LROOTS=',LROOTS
         Write(6,*) ' Nr of roots optimized by super-CI: '//
     &    'NROOTS=',NROOTS
         If (iAll.eq.1) Then
          Write(6,*)' (Equal-weighted)'
         Else
          Write(6,*)' Weights:'
          Do i1=1,NROOTS,10
           i2=min(NROOTS,i1+9)
           write(6,'(1x,10f8.4)')(WEIGHT(i),i=i1,i2)
          End Do
         End If
        End If
       Call ChkIfKey_m()
      ELSE
!AMS - Read the CIROOT information from the old RASSCF run:
       If (DBG) Write(6,*) ' CIROOTS command was not given.'
       If (DBG) Write(6,*) 'I set CIREstart to true.'
       KeyCIRE=.TRUE.
      END IF
*---  Process RLXR command --------------------------------------------*
      If(KeyRLXR) Then
       Call SetPos_m(LUInput,'RLXR',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' O.K. reading after keyword RLXR.'
       Read(LUInput,*,End=9910,Err=9920) IRLXROOT
       ReadStatus=' O.K. reading after keyword RLXR.'
       Call ChkIfKey_m()
       If (DBG) Then
        Write(6,*) ' RLXROOT command was given.'
        Write(6,*) ' State for SLAPAF to handle: ',IRLXROOT
       End If
       If (.not.Any(iRoot(1:LROOTS).eq.IRLXROOT)) Then
        Write(6,*) ' The selected root is not among those available.'
        Call AbEnd()
       End If
      End If
*
*
CIgorS 29-4-2010 Begin
*---  Process MDRL command --------------------------------------------*
      If(KeyMDRL) Then
         If(KeyRLXR) Then
            Write(6,*) ' RLXROOT keyword was given before MDRLXROOT.'
            Write(6,*) ' Since these keywords are mutually exclusive'
            Write(6,*) ' please check the input and read the manual.'
            GoTo 9910
         End If
         Call SetPos_m(LUInput,'MDRL',Line,iRc)
         If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
         Call Qpg_iScalar('Relax CASSCF root',RlxRCheck)
         If (RlxRCheck) THEN
            CALL Get_iScalar('Relax CASSCF root',iRlxRoot)
            If (DBG) Then
               Write(6,*) ' An existing relax root was found.'
               Write(6,*) ' The MDRLXROOT value is ignored.'
            End If
         Else
            ReadStatus=' Failure reading relaxroot number IRLXROOT.'
            Read(LUInput,*,End=9910,Err=9920) IRLXROOT
            ReadStatus=' O.K. after reading relaxroot number IRLXROOT.'
            Call ChkIfKey_m()
            If (.not.Any(iRoot(1:LROOTS).eq.IRLXROOT)) Then
             Write(6,*) ' The selected root is not among those'//
     &                  ' available.'
             Call AbEnd()
            End If
         End If
         If (DBG) Then
            Write(6,*) ' MDRLxroot command was given.'
            Write(6,*) ' DYNAMIX will follow the root: ',IRLXROOT
         End If
      End If
CIgorS End
*
*---  Process CISE command --- (changed!) -----------------------------*
      If (DBG) Write(6,*) ' Check for CISELECT command.'
      If(KeyCISE) Then
       If (DBG) Then
        Write(6,*) ' CISELECT keyword was given.'
        Write(6,*) ' This input is awkward. Let''s find up'
        Write(6,*) ' a better way to do things.'
       End If
       ICICH=1
       Call SetPos_m(LUInput,'CISE',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Do i=1,NROOTS
         ReadStatus=' Failure reading after CISELECT keyword.'
         Read(LUInput,*,End=9910,Err=9920) kRef
         ReadStatus=' O.K. reading after CISELECT keyword.'
         If(kRef.gt.mxRef) Then
            Call WarningMessage(1,'CISElect input is wrong.')
            Write(LF,*) 'Number of CSF''s in CiSelect is out of bounds'
            Write(LF,'(a,i3,a,i3)') 'Specified:',kRef,', Max is',mxRef
            Write(LF,'(a,i3)') 'Standard fixup, value set to',mxRef
            kRef=mxRef
         End If
         ReadStatus=' Failure reading after CISELECT keyword.'
         Read(LUInput,*,End=9910,Err=9920) (ICI(i,iRef),iRef=1,kRef)
         Read(LUInput,*,End=9910,Err=9920) (CCI(i,iRef),iRef=1,kRef)
         ReadStatus=' O.K. reading after CISELECT keyword.'
         dSum=0.0d0
         Do iRef=1,kRef
            dSum=dSum+CCI(i,iRef)**2
         End Do
         Do iRef=1,kRef
            CCI(i,iRef)=CCI(i,iRef)/sqrt(dSum)
         End Do
         Do iRef=kRef+1,mxRef
            CCI(i,iRef)=0.0d0
            ICI(i,iRef)=0
         End Do
       End Do
       Call ChkIfKey_m()
      End If
*
* =========   Input source for orbitals: =============================*
* INVEC=0 is used to indicate if any source of orbitals has been
* identified.
* iOrbData=0 is used to indicate if any source of orbital type
* information -- inactive, ras1, etc -- has been identified.


      If (DBG) Write(6,*) ' Check for input source of orbitals.'
      If (DBG) Then
       Write(6,*)'KeyCORE,KeyJOBI,KeyLUMO:',KeyCORE,KeyJOBI,KeyLUMO
      End If
* Handle multiple keywords:
      If (KeyCORE) Then
        KeyJOBI=.false.
        KeyLUMO=.false.
      Else If (KeyJOBI) Then
        KeyLUMO=.false.
      End if
* Only one of these should have been selected.
      If (DBG) Then
       Write(6,*)'KeyCORE,KeyJOBI,KeyLUMO:',KeyCORE,KeyJOBI,KeyLUMO
      End If

* CORE is probably becoming obsolete.
*---  Process CORE command --------------------------------------------*
      Continue
      If (KeyCORE) Then
       If (DBG) Write(6,*)' CORE command was used.'
        IF (IPRLEV.ge.VERBOSE) Write(LF,*)
     &     ' Start orbitals will be computed from scratch.'
       INVEC=1
      End If

*---  Process JOBI command --------------------------------------------*
      If (KeyJOBI) Then
       If (DBG) Write(6,*)' JOBIPH command was used.'
       Call f_Inquire('JOBOLD',lExists)
       If (lExists) Then
        IF (IPRLEV.ge.VERBOSE) Write(LF,*)
     &     ' Orbitals will be taken from the old jobiph file.'
        INVEC=3
       Else
        Call f_Inquire(IPHNAME,lExists)
        If (lExists) Then
         IF (IPRLEV.ge.VERBOSE) Write(LF,*)
     &     ' Orbitals will be taken from the jobiph file.'
         INVEC=3
        End If
       End If
       If (INVEC.eq.0) Then
        If (IPRLEV.ge.TERSE) Then
          Call WarningMessage(2,'JOBIPH input is wrong.')
         Write(6,*)' Keyword JOBIPH was used, but the ''JOBOLD'' file'
          Write(6,*)' does not exist. The ''JOBIPH'' file named'
          Write(6,*) IPHNAME
          Write(6,*)'also does not exist. This is a fatal error.'
          GoTo 9930
        End If
       End If
      End If

*---  Process H5OR command --------------------------------------------*
      If (KeyH5OR) Then
#ifdef _HDF5_
        KeyLUMO=.false.
        KeyTYPE=.false.
        iOverwr=0
        mh5id = mh5_open_file_r(StartOrbFile)
*     read basic attributes
        call mh5_fetch_attr(mh5id, 'MOLCAS_MODULE', prgm)
        call mh5_fetch_attr(mh5id, 'NSYM', NSYM_L)
        if (nsym.ne.nsym_l) then
          write (6,*) 'Number of symmetries on HDF5 file does not'
          write (6,*) 'match the number of symmetries on the'
          write (6,*) 'RunFile, calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if
        call mh5_fetch_attr(mh5id, 'NBAS', NBAS_L)
        ierr=0
        do isym=1,nsym
          if (nbas(isym).ne.nbas_l(isym)) ierr=1
        end do
        if (ierr.eq.1) then
          write (6,*) 'Number of basis functions on HDF5 file does not'
          write (6,*) 'match the number of basis functions on the'
          write (6,*) 'RunFile, calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if
*     orbitals available?
        if (mh5_exists_dset(mh5id, 'MOCOEF')) then
          inVec=4
        end if
*     typeindex data available?
        if (mh5_exists_dset(mh5id, 'TYPEINDEX')) then
          iOrbData=3
          call mma_allocate(typestring, sum(nbas(1:nsym)))
          call mh5_fetch_dset(mh5id, 'TYPEINDEX', typestring)
          iStart=1
          Do iSym=1,nSym_L
            call tpstr2orb(typestring(iStart),nbas_l(isym),
     $              nFro_L(iSym),nISh_L(iSym),
     $              NRS1_L(iSym),NRS2_L(iSym),NRS3_L(iSym),
     $              nSSh_L(iSym),nDel_L(iSym))
            iStart=iStart+nBas_L(iSym)
          End Do
          call mma_deallocate(typestring)
        end if
*     CI requested?
        if (KeyH5CI) then
          if (         mh5_exists_dset(mh5id, 'CICOEF')
     $            .and.mh5_exists_attr(mh5id, 'LROOTS')) then
            call mh5_fetch_attr(mh5id, 'LROOTS', lroots_l)
            if (lroots_l.ne.lroots) then
              write (6,*) 'Number of CI roots on file does not'
              write (6,*) 'match the number requested by the user,'
              write (6,*) 'so no CI vectors will be read from'
              write (6,*) StartOrbFile
              iCIRST=0
            else
              iCIRST=1
            end if
          else
            write (6,*) 'The required fields CICOEF/LROOTS are'
            write (6,*) 'missing from the HDF5 file supplied by'
            write (6,*) 'the user. As a result, to continue,'
            write (6,*) 'no CI vectors will be read from'
            write (6,*) StartOrbFile
            iCIRST=0
          end if
        end if
        call mh5_close_file(mh5id)
#else
        write (6,*) 'The format of the start orbital file was'
        write (6,*) 'specified by the user as HDF5, but this'
        write (6,*) 'is not implemented in this installation.'
        call Quit(_RC_INPUT_ERROR_)
#endif
      End If

*---  Process LUMO command --------------------------------------------*
      InfoLbl='        '
      If (KeyLUMO) Then
       If (DBG) Then
          Write(6,*)' LUMORB command was used.'
          Write(6,*)' Name of orbital file, StartOrbFile='
          Write(6,*) StartOrbFile
       End If
       Call ChkVec(StartOrbFile,inporb_version,NSYM_L,
     &                                   NBAS_L,NORB_L,InfoLbl,iRc1)
       If(iRc1.ne._RC_ALL_IS_WELL_) Then
        If (IPRLEV.ge.TERSE) Then
         Call WarningMessage(1,'LUMORB input error.')
         Write(6,*)' Keyword LUMORB used with file name StartOrbFile='
         Write(6,*) StartOrbFile
         Write(6,*)' but that file cannot be used. Perhaps it does'
         Write(6,*)' not exist?'
        End If
        iOverWr=1
       Else
        If (DBG) Then
         Write(6,*)' The file may be used to read input orbitals.'
         Write(6,*)' It is of type INPORB version ',inporb_version
         Write(6,*)' The information label is: ',InfoLbl
        End If
* Check that symmetry and basis function information match the runfile:
        IERR=0
        If(NSYM.ne.NSYM_L) THEN
         IERR=1
        Else
         DO ISYM=1,NSYM
          IF(NBAS(ISYM).ne.NBAS_L(ISYM))IERR=2
         END DO
        End If
        IF (IERR.ne. 0) Then
         Call WarningMessage(2,'Unusable start orbital file.')
         Write(6,*)' ERROR: Start orbital file name is '
         Write(6,*) StartOrbFile
         Write(6,*)' That file is a valid orbital file.'
         Write(6,*)' Version:',inporb_version
         IF(IERR.eq.1) Then
          Write(6,*)' Nr of symmetries NSYM_L=',NSYM_L
          Write(6,*)' However, the runfile claims that number'
          Write(6,*)' of symmetries is NSYM=',NSYM
          Write(6,*)' This mismatch is a serious error.'
          Write(6,*)' Calculation must stop -- Sorry!'
          GOTO 9930
         Else
          Write(6,*)' Start orbital file says nr of basis functions'
          Write(6,*)' in each symmetry is:'
          Write(6,'(1x,8I5)')(NBAS_L(I),I=1,NSYM)
          Write(6,*)' but RUNFILE claims it is'
          Write(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
          Write(6,*)' This mismatch is a serious error.'
          Write(6,*)' Calculation must stop -- Sorry!'
          GOTO 9930
         End If
        End If

* This also implies that information on orbital types could be
* taken from typeindex on orbital file:
        If( index(InfoLbl,'i').gt.0  .or. index(InfoLbl,'I').gt.0) Then
          iOrbData=3
          iOverWr=0
          If (DBG) Then
           Write(6,*)' This means we may take orbital specifications'
           Write(6,*)' from the file, so set iOrbData=3, iOverWr=0.'
           Write(6,*)' The orbital spaces are read from typeindex.'
          End If
* We will also take the opportunity to find the orbital spaces size
* according to typeindex, for possible need below:'
          Call GetMem('TypeIdx','Allo','Inte',ipType,mxOrb)
          LuStartOrb=19
          Call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,
     &            Dummy,Dummy,Dummy,iWork(ipType),myTitle,0,iErr)
          call tpidx2orb(NSYM_L,NBAS_L,
     $            iWork(ipType),
     $            NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
          Call GetMem('TypeIdx','Free','Inte',ipType,mxOrb)
        Else
          iOverWr=1
        End If
       End If
       INVEC=2
      End If
      If (DBG) Write(6,*)' The INVEC    code is now',INVEC
      If (DBG) Write(6,*)' The iOrbData code is now',iOrbData
*
* =====================================================================
*
*---  Process TYPEindex command ---------------------------------------*
      If (DBG) Write(6,*)' Was TYPEINDEX requested?'
      If (KeyTYPE) Then
       Call SetPos_m(LUInput,'TYPE',Line,iRc)
       Call ChkIfKey_m()
       If (DBG) Then
        Write(6,*)' TYPEINDEX command was used.'
        Write(6,*)' The size of orbital spaces should be read from'
        Write(6,*)' typeindex in starting orbital file.'
       End If
       If( index(InfoLbl,'i').gt.0  .or. index(InfoLbl,'I').gt.0) Then
         iOrbData=3
         iOverwr=0
         IF(IPRLEV.ge.VERBOSE)
     &    Write(LF,*)' Orbital specification will be taken '//
     &               'from orbital file'
          Call GetMem('TypeIdx','Allo','Inte',ipType,mxOrb)
          LuStartOrb=19
          Call RdVec(StartOrbFile,LuStartOrb,'IA',NSYM_L,NBAS_L,NBAS_L,
     &            Dummy,Dummy,Dummy,iWork(ipType),myTitle,0,iErr)
          call tpidx2orb(NSYM_L,NBAS_L,
     $            iWork(ipType),
     $            NFRO_L,NISH_L,NRS1_L,NRS2_L,NRS3_L,NSSH_L,NDEL_L)
          Call GetMem('TypeIdx','Free','Inte',ipType,mxOrb)
         IERR=0
         IF (NSYM_L.ne.NSYM) IERR=1
         IF(IERR.eq.0) THEN
          DO ISYM=1,NSYM
           IF(NBAS(ISYM).ne.NBAS_L(ISYM))IERR=1
          END DO
         END IF
         IF (IERR.ne.0 .or. DBG) THEN
          Write(LF,*)
          Write(LF,'(6X,A)')'Specifications read from runfile:'
          Write(LF,'(6X,A)')'----------------------------------------'
          Write(LF,*)
          Write(LF,'(6X,A,T47,8I4)') 'Symmetry species',
     &                                (iSym,iSym=1,NSYM)
          Write(LF,'(6X,A,T47,8I4)') 'Number of basis functions',
     &                                (NBAS(iSym),iSym=1,NSYM)
          Write(LF,*)
          Write(LF,'(6X,A)')'Specifications read from orbital file:'
          Write(LF,'(6X,A)')'----------------------------------------'
          Write(LF,*)
          Write(LF,'(6X,A,T47,8I4)') 'Symmetry species',
     &                                (iSym,iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'Frozen orbitals',
     &                                (NFRO_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'Inactive orbitals',
     &                                (NISH_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'RAS1 orbitals',
     &                                (NRS1_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'RAS2 orbitals',
     &                                (NRS2_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'RAS3 orbitals',
     &                                (NRS3_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'Secondary orbitals',
     &                                (NSSH_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'Deleted orbitals',
     &                                (NDEL_L(iSym),iSym=1,NSYM_L)
          Write(LF,'(6X,A,T47,8I4)') 'Number of basis functions',
     &                                (NBAS_L(iSym),iSym=1,NSYM_L)
          Write(LF,*)
         END IF
         IF (IERR.ne.0 .and.IPRLEV.ge.TERSE) THEN
          Write(6,*)' Orbital specifications were to be read from'
          Write(6,*)' orbital file, but there is mismatch with'
          Write(6,*)' some data on the runfile!'
          Write(6,*)' Orbital file name is:'
          Write(6,*) StartOrbFile
          iOrbData=0
         ELSE
          DO ISYM=1,NSYM
           NFRO(ISYM)=NFRO_L(ISYM)
           NISH(ISYM)=NISH_L(ISYM)
           NRS1(ISYM)=NRS1_L(ISYM)
           NRS2(ISYM)=NRS2_L(ISYM)
           NRS3(ISYM)=NRS3_L(ISYM)
           NSSH(ISYM)=NSSH_L(ISYM)
           NDEL(ISYM)=NDEL_L(ISYM)
          END DO
         END IF
       Else
        If (DBG) Write(6,*)' Typeindex cannot be read!'
       End If
      End If
      If (DBG) Write(6,*)' The iOrbData code is now',iOrbData

* =======================================================================
* Explicit orbital sizes input:
* Save a copy of current iorbdata first:
      iod_save=iorbdata
*---  Process FROZ command --------------------------------------------*
      If (KeyFROZ) Then
       If (DBG) Write(6,*) ' FROZEN keyword was given.'
       Call SetPos_m(LUInput,'FROZ',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after FROZEN keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NFRO(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. reading after FROZEN keyword.'
       Call Get_iScalar('nSym',i)
       Call Put_iArray('nFro',nFro,i)
       If (DBG) Then
        Write(6,*) ' Nr of Frozen orbitals requested:'
        Write(6,'(1x,8i5)')(NFRO(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If
*---  Process INAC command --------------------------------------------*
      If (KeyINAC) Then
       If (DBG) Write(6,*) ' INACTIVE keyword was given.'
       Call SetPos_m(LUInput,'INAC',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after INACTIVE keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NISH(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. reading after INACTIVE keyword.'
       If (DBG) Then
        Write(6,*) ' Nr of Inactive orbitals requested:'
        Write(6,'(1x,8i5)')(NISH(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If
*
*---  Process RAS1 command --------------------------------------------*
      If (KeyRAS1) Then
       If (DBG) Write(6,*) ' RAS1 keyword was given.'
       Call SetPos_m(LUInput,'RAS1',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after RAS1 keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NRS1(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. reading after RAS1 keyword.'
       If (DBG) Then
        Write(6,*) ' Nr of RAS1 orbitals requested:'
        Write(6,'(1x,8i5)')(NRS1(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If
*
*---  Process RAS2 command --------------------------------------------*
      If (KeyRAS2) Then
       If (DBG) Write(6,*) ' RAS2 keyword was given.'
       Call SetPos_m(LUInput,'RAS2',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after RAS2 keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NRS2(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. reading after RAS2 keyword.'
       If (DBG) Then
        Write(6,*) ' Nr of Ras2 orbitals requested:'
        Write(6,'(1x,8i5)')(NRS2(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If
*
*---  Process RAS3 command --------------------------------------------*
      If (KeyRAS3) Then
       If (DBG) Write(6,*) ' RAS3 keyword was given.'
       Call SetPos_m(LUInput,'RAS3',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after RAS3 keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NRS3(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. reading after RAS3 keyword.'
       If (DBG) Then
        Write(6,*) ' Nr of Ras3 orbitals requested:'
        Write(6,'(1x,8i5)')(NRS3(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If
*
*---  Process GASS command --------------------------------------------*
      If (KeyGASS) Then
       if(DBG) write(6,*) 'GAS is actived'
       call setpos_m(luinput,'GASS',line,irc)
       if(irc.ne._RC_ALL_IS_WELL_) goto 9810
       read(luinput,*,end=9910,err=9920) NGAS
       do igas=1,ngas
         read(luinput,*,end=9910,err=9920)
     &    (ngssh(igas,isym),isym=1,nsym)
         read(luinput,*,end=9910,err=9920)
     &    (igsoccx(igas,mm),mm=1,2)
       end do
       iDoGas = .true.
       iorbdata=1
      end if
*
*---  Process DELE command --------------------------------------------*
      If (KeyDELE) Then
       If (DBG) Write(6,*) ' DELETED keyword was given.'
       Call SetPos_m(LUInput,'DELE',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading after DELETED keyword.'
       Read(LUInput,*,End=9910,Err=9920) (NDEL(ISYM),ISYM=1,NSYM)
       ReadStatus=' O.K. after reading after DELETED keyword.'
       If (DBG) Then
        Write(6,*) ' Nr of Deleted orbitals requested:'
        Write(6,'(1x,8i5)')(NDEL(i),i=1,NSYM)
       End If
       IORBDATA=1
      End If

      IF (IORBDATA.eq.1 .and. IPRLEV.ge.VERBOSE) THEN
       Write(LF,*)' Orbital specification was read from input.'
      END IF
      IF (IOD_SAVE.eq.3 .and. IORBDATA.eq.1) THEN
* See if the input matches the values on file
       IERR=0
       DO ISYM=1,NSYM
        IF(NFRO(ISYM).ne.NFRO_L(ISYM))IERR=1
        IF(NISH(ISYM).ne.NISH_L(ISYM))IERR=1
        IF(NRS1(ISYM).ne.NRS1_L(ISYM))IERR=1
        IF(NRS2(ISYM).ne.NRS2_L(ISYM))IERR=1
        IF(NRS3(ISYM).ne.NRS3_L(ISYM))IERR=1
        IF(NSSH(ISYM).ne.NSSH_L(ISYM))IERR=1
        IF(NDEL(ISYM).ne.NDEL_L(ISYM))IERR=1
       END DO
       IF(IERR.eq.0) THEN
        IF (IORBDATA.eq.1 .and. IPRLEV.ge.VERBOSE) THEN
         Write(LF,*)' However, input matches the typeindex on the'
         Write(LF,*)' starting orbitals file. Therefore, accept'
         Write(LF,*)' the typeindex information for sorting.'
         iOrbData=3
         iOverWr=0
        END IF
       END IF
      END IF
      IF(IORBDATA.eq.3) THEN
       DO ISYM=1,NSYM
        NFRO(ISYM)=NFRO_L(ISYM)
        NISH(ISYM)=NISH_L(ISYM)
        NRS1(ISYM)=NRS1_L(ISYM)
        NRS2(ISYM)=NRS2_L(ISYM)
        NRS3(ISYM)=NRS3_L(ISYM)
        NSSH(ISYM)=NSSH_L(ISYM)
        NDEL(ISYM)=NDEL_L(ISYM)
       END DO
      END IF
* =======================================================================
* If IORBDATA is still 0, lets hope there is information on the runfile.
* Exception: If this is a CIRESTART, it must be taken from the JOBIPH
* (or JOBOLD) file.
      iprlev=insane
      IF (IORBDATA.EQ.0) THEN
       IF(IPRLEV.ge.VERBOSE)
     &          Write(LF,*)' No explicit orbital specs in user input.'
       IF(KeyCIRE) Then
!        IF(IPRLEV.ge.VERBOSE) Then
!         Write(LF,*)' This is a CIRESTART case, so take them from'
!         Write(LF,*)' the JOBIPH or JOBOLD file.'
!        End If
!        IORBDATA=2
!        IF(IPRLEV.ge.VERBOSE)
!     &    Write(LF,*)' Orbital specs taken from JOBIPH or JOBOLD.'
        IAD19=0
!        iJOB=0
!        Call f_Inquire('JOBOLD',lExists)
!        If (lExists) iJOB=1
!        if(JOBOLD.le.0) Then
!          JOBOLD=20
!        end if
!        If (iJOB.eq.1) Then
!           Call DaName(JOBOLD,'JOBOLD')
!        Else
!           If (IPRLEV.ge.TERSE) then
!              Call WarningMessage(1,'JOBOLD not found, using JOBIPH.')
!           End If
!           If (JOBIPH.gt.0) Then
!              JOBOLD=JOBIPH
!           Else
!              Call DaName(JOBOLD,'JOBIPH')
!           End If
!        End If
!        Call IDaFile(JOBOLD,2,IADR19,10,IAD19)
        Call IDaFile(JOBIPH,2,IADR19,10,IAD19)
        lll = 1
        lll = MAX(lll,mxSym)
        lll = MAX(lll,mxOrb)
        lll = MAX(lll,RtoI)
        lll = MAX(lll,4*2*mxOrb/ItoB)
        lll = MAX(lll,2*72/ItoB)
        lll = MAX(lll,RtoI*mxRoot)
        CALL GETMEM('JOBOLD','ALLO','INTEGER',lJobH,lll)
* PAM Jan 2014 -- do not take POTNUC from JOBIPH; take it directly
* from runfile, where it was stored by seward.
        iAd19=iAdr19(1)
!        CALL WR_RASSCF_Info(JobOld,2,iAd19,NACTEL,ISPIN,NSYM,LSYM,
        CALL WR_RASSCF_Info(JobIPH,2,iAd19,NACTEL,ISPIN,NSYM,LSYM,
     &                      NFRO,NISH,NASH,NDEL,NBAS,
     &                      mxSym,iWork(lJobH),4*2*mxOrb,NCONF,
     &                      iWork(lJobH),2*72,JobTit,4*18*mxTit,
     &                      POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,
     &                      NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
        CALL GETMEM('JOBOLD','FREE','INTEGER',lJobH,lll)
       ELSE
        IF(IPRLEV.ge.VERBOSE) Then
         Write(LF,*)' This is not a CIRESTART case, so take them from'
         Write(LF,*)' the RUNFILE.'
        End If
        Call qpg_iArray('nAsh',lExists,nItems)
        IF (lExists .and. nItems.eq.nSym) THEN
         Call Get_iArray('nFro',nFro,nSym)
         Call Get_iArray('nISh',nISh,nSym)
         Call Get_iArray('nASh',nRS2,nSym)
         Call Get_iArray('nDel',nDel,nSym)
         IORBDATA=2
        END IF
       END IF
      End If

!AMS - this may be closing either JOBOLD or JOBIPH. Close only JOBOLD.
      IF(JOBOLD>0) then
        IF(JOBOLD.ne.JOBIPH) THEN
            Call DaClos(JOBOLD)
        END IF
      end if


!AMS - make sure we change to a different JOBIPH file - we don't want to
!overwrite any existing JOBIPH file.
!
!Close the old JOBIPH file
      if(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
      end if
!Rename JOBIPH file, and open it.
!      IPHNAME='JOBPDFT'
      JOBIPH=IsFreeUnit(15)
      CALL DANAME(JOBIPH,IPHNAME)

*---  complete orbital specifications ---------------------------------*
      Do iSym=1,nSym
        if(.not.iDoGas)then
          nash(isym)=nrs1(isym)+nrs2(isym)+nrs3(isym)
        else
          NASH(ISYM)=SUM(NGSSH(1:NGAS,ISYM))
        end if
        NORB(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
        NSSH(ISYM)=NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
      End Do
*---  Related data for sizes, etc.
      NTOT=0
      NTOT1=0
      NTOT2=0
      NO2M=0
      NISHT=0
      NASHT=0
      NDELT=0
      NFROT=0
      NSEC=0
      NORBT=0
      NTOT3=0
      NTOTSP=0
      NTOT4=0
      NRS1T=0 ! for RASSCF
      NRS2T=0
      NRS3T=0
c      Call FZero(NGSSH_tot,ngas)
c      do igas=1,ngas
c        NGSSH_tot(igas) = SUM(NGSSH(IGAS,1:NSYM))
c      end do
c      if(dbg) then
c        write(6,*) 'NGSSH_tot(igas):'
c        write(6,*) (NGSSH_tot(igas),igas=1,ngas)
c      end if
      DO ISYM=1,NSYM
         NTOT=NTOT+NBAS(ISYM)
         NTOT1=NTOT1+NBAS(ISYM)*(NBAS(ISYM)+1)/2
         NTOT2=NTOT2+NBAS(ISYM)**2
         NO2M=MAX(NO2M,NBAS(ISYM)**2)
         NRS1T=NRS1T+NRS1(ISYM)  ! for RAS
         NRS2T=NRS2T+NRS2(ISYM)
         NRS3T=NRS3T+NRS3(ISYM)
         NFROT=NFROT+NFRO(ISYM)
         NISHT=NISHT+NISH(ISYM)
         NASHT=NASHT+NASH(ISYM)
         NDELT=NDELT+NDEL(ISYM)
         NSEC=NSEC+NSSH(ISYM)
         NORBT=NORBT+NORB(ISYM)
         NTOT3=NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
         NTOTSP=NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
         NTOT4=NTOT4+NORB(ISYM)**2
      END DO
      NACPAR=(NASHT+NASHT**2)/2
      NACPR2=(NACPAR+NACPAR**2)/2
* NASHT is called NAC in some places:
      NAC=NASHT
* Same, NISHT, NIN:
      NIN=NISHT
      NFR=NFROT
      If (DBG) Write(6,*)' The iOrbData code is now',iOrbData
* =======================================================================
* Compute effective nuclear charge.
* Identical to nr of protons for conventional basis sets only, not ECP.
      Call Get_iScalar('Unique atoms',nNuc)
      Call GetMem('EffNChrg','Allo','Real',ipENC,nNuc)
      Call Get_dArray('Effective nuclear Charge',Work(ipENC),nNuc)
      TEffNChrg=0.0D0
      Call GetMem('nStab','Allo','Inte',ipStab,nNuc)
      Call Get_iArray('nStab',iWork(ipStab),nNuc)
      do i=1,nNuc
       TEffNChrg=TEffNChrg+Work(ipENC-1+i)*DBLE(nSym/iWork(ipStab-1+i))
      end do
      Call GetMem('nStab','Free','Inte',ipStab,nNuc)
      Call GetMem('EffNChrg','Free','Real',ipENC,nNuc)
      If (DBG) Write(6,*)
     &             ' Effective nuclear charge is TEffNChrg=',TEffNChrg
      TotChrg=0.0D0
      If (DBG) Write(6,*)' Set TotChrg=',TotChrg
*---  Process NACT command --------------------------------------------*
      IF(KeyNACT) Then
C Cannot set the number of inactive orbitals if explicitly given or
C or if there is symmetry
       If(KeyCHAR.and.(KeyINAC.or.(NSYM.gt.1))) Then
        If(IPRLEV.ge.TERSE) Then
         Write(6,*)' CHARGE and NACTEL are only compatible if'
         Write(6,*)' INACTIVE is not given and with no symmetry (C1)'
         Write(6,*)' The CHARGE command will be ignored.'
        End If
        KeyCHAR=.false.
       End If
       If (DBG) Write(6,*) ' NACTEL keyword was given.'
       Call SetPos_m(LUInput,'NACT',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Line=Get_Ln(LUInput)
       ReadStatus=' Failure reading data after NACTEL keyword.'
       Read(Line,*,Err=9920,End=9920) NACTEL,NHOLE1,NELEC3
       ReadStatus=' O.K. after reading data after NACTEL keyword.'
      If (DBG) Write(6,*)' NACTEL,NHOLE1,NELEC3:',NACTEL,NHOLE1,NELEC3
C Only set total charge here if not explicitly given
       If(.not.KeyCHAR) Then
        TotChrg=TEffNChrg-DBLE(2.0D0*(NISHT+NFROT)+NACTEL)
        If (DBG) Write(6,*) ' TotChrg=',TotChrg
       EndIf
       Call Put_iScalar('nActel',NACTEL)
       Call ChkIfKey_m()
      End If
*---  Process CHARGE command --------------------------------------------*
      IF(KeyCHAR) Then
       If (DBG) Write(6,*) ' CHARGE command was given.'
       Call SetPos_m(LUInput,'CHAR',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Line=Get_Ln(LUInput)
       ReadStatus=' Failure reading charge after CHARGE keyword.'
       Read(Line,*,Err=9920) NCHRG
       ReadStatus=' O.K. reading charge after CHARGE keyword.'
       If (DBG) Write(6,*) ' Total charge is ',NCHRG
       TotChrg=DBLE(NCHRG)
       Call ChkIfKey_m()
C If both CHAR and NACT where given, set the inactive and secondary
C orbitals accordingly
       If (KeyNACT) Then
         NISHT_old=NISHT
         NISHT=(INT(TEffNChrg-TotChrg+0.5D0)-NACTEL)/2-NFROT
         NISH(1)=NISHT
         NIN=NISHT
         NSEC=NSEC-(NISHT-NISHT_old)
         NSSH(1)=NSEC
       End If
      End If
* The NINT function is unreliable on Cygwin gfortran, use INT:
* Nr of electrons should be positive integer, so this is probably safe:
      If (DBG) Then
        write(6,*)'NACTEL must be either read from JOBIPH/JOBOLD'
        write(6,*)'or set using the NACT keyword.'
        write(6,*) 'NACTEL=',NACTEL
      end if
!      NACTEL=INT(TEffNChrg-TotChrg+0.5D0)-2*(NISHT+NFROT)
!      Call Put_iScalar('nActel',NACTEL)
!      If (DBG) Then
!        write(6,*)' Compute NActEl from  other data:'
!        write(6,*)'     TEffNChrg=',TEffNChrg
!        write(6,*)'       TotChrg=',TotChrg
!        write(6,*)'       NISHT  =',NISHT
!        write(6,*)'       NFROT  =',NFROT
!        write(6,*)' Resulting NActEl=',NActEl
!      End If
*---  Process RASSCF command --------------------------------------------*
      IF(KeyRASS) Then
       If (DBG) Write(6,*) ' RASSCF keyword was given.'
       Call SetPos_m(LUInput,'RASS',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Line=Get_Ln(LUInput)
       ReadStatus=' Failure reading data after RASSCF keyword.'
       Read(Line,*,Err=9920,End=9920) NHOLE1,NELEC3
       ReadStatus=' O.K. reading data after RASSCF keyword.'
       Call ChkIfKey_m()
      End If
* =======================================================================
*---  Process SPIN command --------------------------------------------*
      If (DBG) Write(6,*)' Determine spin value:'
      IF(KeySPIN) Then
       If (DBG) Write(6,*) ' SPIN command was given.'
       Call SetPos_m(LUInput,'SPIN',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Line=Get_Ln(LUInput)
       ReadStatus=' Failure reading spin after SPIN keyword.'
       Read(Line,*,Err=9920) ISPIN
       ReadStatus=' O.K. reading spin after SPIN keyword.'
       If (DBG) Then
         Write(6,*) ' Spin multiplicity is ',ISPIN
         S=0.5D0*DBLE(ISPIN-1)
         Write(6,*) '  i.e. SPIN=',S
       End If
       Call ChkIfKey_m()
* If ISPIN has not been set, use some value that may have been set
* in runfile (e.g., UHF, or previous calculation)
       ELSE IF(KeyCIRE) Then
        If (DBG) Write(6,*) ' SPIN is taken from JOBOLD (CIRestart=T)'
!        Write(6,*)' ISPIN=',ISPIN

       ELSE
        If (DBG) Write(6,*) ' No SPIN command was given.'
        Call qpg_iscalar('ISPIN',lExists)
        if(lExists) then
          Call get_iscalar('ISPIN',ISPIN)
          If (DBG) Write(6,*) ' Runfile has ISPIN=',ISPIN
        else
          If (DBG) Write(6,*) ' Runfile does not know ISPIN.'
          Call qpg_dscalar('UHFSPIN',lExists)
          if(lExists) then
           Call get_dscalar('UHFSPIN',SUHF)
           ISPIN=NINT(1.0D0+2.0D0*SUHF)
           If (DBG) Write(6,*) ' Runfile has UHFSPIN=',SUHF
          else
* or, last chance fallback, guess on singlet or doublet:
           If (DBG) Write(6,*) ' Runfile does not know UHFSPIN.'
           ISPIN=1
           IF(NACTEL.ne.2*(NACTEL/2)) ISPIN=2
           If(IPRLEV.ge.TERSE) Then
            Call WarningMessage(1,'Had to guess the spin.')
            Write(6,*)' Warning: no input and no reliable source'
            Write(6,*)' for the spin multiplicity.'
            Write(6,*)' Guess ISPIN=',ISPIN
           End If
          end if
        end if
      END IF
      Call put_iscalar('ISPIN',ISPIN)
* If spin is zero, do not compute and print spin density:
      IF (ISPIN.EQ.1) THEN
        ISPDEN=0
      END IF
* =======================================================================
      IF(KeySYMM) Then
*---  Process SYMM command --------------------------------------------*
       If (DBG) Write(6,*) ' SYMMETRY command was given.'
       Call SetPos_m(LUInput,'SYMM',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Line=Get_Ln(LUInput)
       ReadStatus=' Failure reading symmetry index after SYMM keyword.'
       Read(Line,*,Err=9920) LSYM
       ReadStatus=' O.K. reading symmetry index after SYMM keyword.'
       If (DBG) Write(6,*) ' State symmetry index ',LSYM
       Call ChkIfKey_m()
* If LSYM has not been set, normally it should be defaulted to 1.
* Exception: if this is a high-spin OS case, these often require LSYM.ne.1:
       ELSE IF(KeyCIRE) Then
        If (DBG) Write(6,*) ' SYMM is taken from JOBOLD (CIRestart=T)'
       ELSE
        LSYM=1
        IF(ISPIN.eq.NASHT+1) THEN
         DO ISYM=1,NSYM
          NA=NASH(ISYM)
          IF(NA.ne.2*(NA/2)) LSYM=MUL(LSYM,ISYM)
         END DO
        END IF
      END IF
      Call put_iscalar('LSYM',LSYM)
      If (DBG) Write(6,*)' State symmetry LSYM=',LSYM
*
* =======================================================================
*---  Process CIRE command --------------------------------------------*
      If (KeyCIRE) Then
       If (DBG) Write(6,*) ' CIRESTART keyword was given.'
       ICIRST=1
      End If
*
*---  Process HOME command (root homing in SXCI part)--------------------*
      If (KeyHOME) Then
       SXSEL='HOMING  '
       If (DBG) Write(6,*) ' HOME (Root Homing) keyword was given.'
        Call SetPos_m(LUInput,'HOME',Line,iRc)
        Call ChkIfKey_m()
      End If

*---  Process SUPS command ---
      If (KeySUPS) Then
       If (DBG) Write(6,*) ' SUPS (Supersymmetry) keyword was given.'
       Call SetPos_m(LUInput,'SUPS',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Call GetMem('Temp1','Allo','Inte',ipTemp1,mxOrb)
       ISUPSM=1
       iOffset=0
       Do iSym=1,nSym
          ReadStatus=' Failure reading data following SUPS keyword.'
          Read(LUInput,*,End=9910,Err=9920) nGrp
          ReadStatus=' O.K. after reading data following SUPS keyword.'
          Do iGrp=1,nGrp
             Call RdSups_m(LUInput,kOrb,iWork(ipTemp1))
             Do iOrb=1,kOrb
                IXSYM(iWork(ipTemp1+iOrb-1)+iOffset)=iGrp
             End Do
          End Do
          iOffset=iOffset+nBas(iSym)
       End Do
       Call GetMem('Temp1','Free','Inte',ipTemp1,mxOrb)
* (SVC) If both ALTER and SUPS keyword has been used, then change the IXSYM
* arrays according to the changed orbital ordering given in ALTER.
       Do iAlter=1,NAlter
         iChng1=IXSYM(iMAlter(iAlter,1))
         iChng2=IXSYM(iMAlter(iAlter,2))
         IXSYM(iMAlter(iAlter,1))=iChng2
         IXSYM(iMAlter(iAlter,2))=iChng1
       End Do
       Call ChkIfKey_m()
      End If
*
*
*---  Process CLEA command ---
      Continue
      If (KeyCLEA) Then
       If (DBG) Write(6,*) ' CLEAN (Orbital Cleaning) keyword.'
       If (DBG) Write(6,*) ' (Awkward input -- replace??).'
       Call SetPos_m(LUInput,'CLEA',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Call GetMem('Temp1','Allo','Inte',ipTemp1,mxOrb)
       Call GetMem('Temp2','Allo','Inte',ipTemp2,mxOrb)
       Call GetMem('Temp3','Allo','Inte',ipTemp3,mxOrb)
       ICLEAN=1
       nClean=0
       Do iSym = 1, nSym
          nClean=nClean+nBas(iSym)**2
       End Do
       Call GetMem('CleanMask','Allo','INTE',ipCleanmask,nClean)
       iOffset = ipCleanMask-1
       Do iSym=1,nSym
         mBas = nBas(iSym)
         Do i = 1,mBas
           ii = (i-1)*mBas
           Do j = 1,mBas
             ij = j+ii+iOffset
             iWork(ij) = 0
           End Do
         End Do
         ReadStatus=' Failure reading data following CLEAN keyword.'
         Read(LUInput,*,End=9910,Err=9920) nGrp
         ReadStatus=' O.K. after reading data following CLEAN keyword.'
         Do iGrp = 1,nGrp
          Call RdSups_m(LUInput,mOrb,iWork(ipTemp1))
          Call RdSups_m(LUInput,nCof,iWork(ipTemp2))
          Call RdSups_m(LUInput,mCof,iWork(ipTemp3))
          Do i = 1,mBas
            ii = (i-1)*mBas
            is_in_Group = 0
            Do j = 1,mOrb
              If ( iWork(ipTemp1+j-1).eq.i ) is_in_Group = 1
            End Do
            If ( is_in_Group.eq.1 ) then
              Do k = 1,nCof
                ij = iWork(ipTemp2+k-1)+ii+iOffset
                iWork(ij) = 1
              End Do
            Else
              Do k = 1,mCof
                ij = iWork(ipTemp3+k-1)+ii+iOffset
                iWork(ij) = 1
              End Do
            End If
          End Do
         End Do
         iOffset = iOffset+mBas*mBas
       End Do
       Call GetMem('Temp1','Free','Inte',ipTemp1,mxOrb)
       Call GetMem('Temp2','Free','Inte',ipTemp2,mxOrb)
       Call GetMem('Temp3','Free','Inte',ipTemp3,mxOrb)
       Call ChkIfKey_m()
      End If

*---  Process CHOL command (Cholesky Default Input, F.Aquilante Sept 04)
      If (KeyCHOL) Then
       If (DBG) Write(6,*) ' CHOLESKY keyword was given.'
       DoCholesky=.True.
       Call SetPos_m(LUInput,'CHOL',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Call Cho_rasscf_rdInp_m(.True.,LUInput)
      End If
*
*---  Process CHOI command (Cholesky Custom Input, F.Aquilante)
* Cholesky with user-defined settings.
      If (KeyCHOI) Then
       If (DBG) Write(6,*) ' CHOINPUT keyword was given.'
       DoCholesky=.True.
       Call SetPos_m(LUInput,'CHOI',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Call Cho_rasscf_rdInp_m(.False.,LUInput)
      End If
*
*---  Process ITER command --------------------------------------------*
      If (KeyITER) Then
       If (DBG) Write(6,*) ' ITERATIONS keyword was given.'
       Call SetPos_m(LUInput,'ITER',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data following ITER keyword.'
       Read(LUInput,*,End=9910,Err=9920) MAXIT,ITMAX
       ReadStatus=' O.K. reading data following ITER keyword.'
       If (DBG) Then
        Write(6,*) ' Max nr of RASSCF (macro) iterations MAXIT=',MAXIT
        Write(6,*) ' Max nr of orbital optimization iter ITMAX=',ITMAX
        If (KeyCION) Then
         Write(6,*)
     &  ' (Irrelevant, since the CIONLY keyword was also given)'
        End If
       End If
       Call ChkIfKey_m()
      End If
*
*---  Process LEVS command --------------------------------------------*
      If (KeyLEVS) Then
       If (DBG) Write(6,*) ' LEVS (Level Shift) keyword was used.'
       Call SetPos_m(LUInput,'LEVS',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading level shift after LEVSHFT keyword.'
       Read(LUInput,*,End=9910,Err=9920) LVSHFT
       ReadStatus=' O.K. reading level shift after LEVSHFT keyword.'
       If (DBG) Write(6,*) ' Level shift LVSHFT (Re*8!!) =',LVSHFT
       Call ChkIfKey_m()
      End If
*
*---  Process THRS command --------------------------------------------*
      If (KeyTHRS) Then
       If (DBG) Write(6,*) ' THRS (Thresholds) command was used.'
       Call SetPos_m(LUInput,'THRS',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading thresholds after THRS keyword.'
       Read(LUInput,*,End=9910,Err=9920) THRE,THRTE,THRSX
       ReadStatus=' O.K. after reading thresholds after THRS keyword.'
*tbp, may 2013: no altering of thre with Cholesky
*tbp   If (DoCholesky.and.IPRLEV.GE.TERSE) then
*tbp     write(6,*)'*** Detected Cholesky or RI/DF calculation'
*tbp     write(6,*)' BUT user specified value of ThrE will be used.'//
*tbp &             ' ThrE= ',THRE
*tbp   End If
       If (DBG) Then
        Write(6,*) ' Threshold for energy change,    THRE =',THRE
        Write(6,*) ' Threshold for orbital rotation, THRTE=',THRTE
        Write(6,*) ' Threshold for max BLB element,  THRSX=',THRSX
       End If
       Call ChkIfKey_m()
      End If
*
*---  Process TIGH command --------------------------------------------*
      If (KeyTIGH) Then
       If (DBG) Write(6,*) ' TIGHT (Tight CI convergence)  used.'
       KTIGHT=1
       Call SetPos_m(LUInput,'TIGH',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after TIGHT keyword.'
       Read(LUInput,*,End=9910,Err=9920) THREN,THFACT
       ReadStatus=' O.K. reading data after TIGHT keyword.'
       If (DBG) Then
        Write(6,*) ' CI energy threshold in 1st iter, THREN=',THREN
        Write(6,*) ' CI Threshold/Energy change =    THFACT=',THFACT
       End If
       Call ChkIfKey_m()
      End If
*
* Use of quasi-newton ci/orbital coupling?
* Commands QUNE or NOQUNE:
      If (KeyNOQU) Then
        NQUNE=0
        If (DBG) Write(6,*) ' NOQUNE keyword: QUNE is disabled.'
        Call SetPos_m(LUInput,'NOQU',Line,iRc)
        Call ChkIfKey_m()
      Else If (KeyQUNE) Then
        If (DBG) Write(6,*) ' QUNE keyword: QUNE is enabled.'
        NQUNE=1
        Call SetPos_m(LUInput,'QUNE',Line,iRc)
        Call ChkIfKey_m()
      Else
* Default is to use QUNE, unless this is some kind of DFT:
       If (KeyKSDF) Then
         NQUNE=0
         If (DBG) Write(6,*) ' DFT calculation: QUNE is disabled.'
       Else
         NQUNE=1
         If (DBG) Write(6,*) ' QUNE is enabled by default.'
       End If
      End If
*
*---  Process AVER command --------------------------------------------*
      If (KeyAVER) Then
       If (DBG) Write(6,*) ' AVER (Symmetry averaging) is used.'
       If (DBG) Write(6,*) ' This is probably becoming obsolete.'
       Call SetPos_m(LUInput,'AVER',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after KAVER keyword.'
       Read(LUInput,*,End=9910,Err=9920) KAVER,(KSYM(I),I=1,2*KAVER)
       ReadStatus=' O.K. after reading data after KAVER keyword.'
       If (KAVER.ne.1 .and. KAVER.ne.2) Then
        If (IPRLEV.GE.TERSE) Then
         write(6,*)' AVER keyword used with inappropriate numbers'
         write(6,*)' and is ignored.'
        End If
        KAVER=0
        KeyAVER=.FALSE.
       End if
       If (DBG) Then
        Write(6,*) ' Averaging (blindly) density matrices of symmetry'
        Write(6,*) ' labels ',KSYM(1),' and ',KSYM(2)
        If (KAVER.eq.2) Then
         Write(6,*) ' and also labels ',KSYM(3),' and ',KSYM(4)
        End If
       End If
       Call ChkIfKey_m()
      End If
*
*---  Process CIMX command --------------------------------------------*
      If (KeyCIMX) Then
       If (DBG) Write(6,*) ' Keyword CIMX was used.'
       Call SetPos_m(LUInput,'CIMX',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after CIMX keyword.'
       Read(LUInput,*,End=9910,Err=9920) MAXJT
       ReadStatus=' O.K. after reading data after CIMX keyword.'
       If (DBG) Write(6,*) ' Max nr of CI iterations MAXJT=',MAXJT
       Call ChkIfKey_m()
      End If
*
*---  Process SDAV command --------------------------------------------*
      If (KeySDAV) Then
       If (DBG) Then
         Write(6,*) ' SDAV (Size of explicit Hamiltonian matrix)'
       End If
       Call SetPos_m(LUInput,'SDAV',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after SDAV keyword.'
       Read(LUInput,*,End=9910,Err=9920) NSEL
       ReadStatus=' O.K. after reading data after SDAV keyword.'
       If (DBG) Write(6,*) ' Size is NSEL=',NSEL
       Call ChkIfKey_m()
      End If
*
*---  Process OFEM commands for Orbital-Free embedding -------------*
*
      If (KeyOFEM) Then
       If (DBG) Then
         Write(6,*) ' OFEM (Orbital-Free Embedding activated)'
       End If
       Do_OFemb=.true.
       Call SetPos_m(LuInput,'OFEM',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after OFEM keyword.'
       Read(LUInput,'(A)',End=9910,Err=9920) OFE_KSDFT
       ReadStatus=' O.K. after reading data after OFEM keyword.'
       Call UpCase(OFE_KSDFT)
       Call LeftAd(OFE_KSDFT)
       write(6,*)
       write(6,*)  '  --------------------------------------'
       write(6,*)  '   Orbital-Free Embedding Calculation'
       write(6,*)  '  --------------------------------------'
       If (OFE_KSDFT(1:4).eq.'LDTF') Then
        write(6,*) '    T_nad potential   : Thomas-Fermi    '
       Else
        write(6,*) '    T_nad potential   : ',OFE_KSDFT(1:4)
       EndIf
       If (KEonly) Then
        write(6,*) '    Exc_nad potential :  None           '
       Else
        write(6,*) '    Exc_nad potential : ',OFE_KSDFT(6:10)
       EndIf
       write(6,*)  '  --------------------------------------'
       write(6,*)
       DFTFOCK='SCF '
      EndIf
      If (KeyFTHA) Then
       Call SetPos_m(LuInput,'FTHA',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after FTHA keyword.'
       Read(LUInput,*,End=9910,Err=9920) ThrFThaw
       ReadStatus=' O.K. after reading data after FTHA keyword.'
      EndIf
      If (KeyDFMD) Then
       Call SetPos_m(LuInput,'DFMD',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after DFMD keyword.'
       Read(LUInput,*,End=9910,Err=9920) dFMD, Xsigma
       ReadStatus=' O.K. after reading data after DFMD keyword.'
c       Write(6,'(A,F6.3)') '  Fract. correl. potent. DFMD=',dFMD
c       write(6,*)          '  --------------------------------------'
       If (dFMD+Xsigma.lt.0.0d0) Then
      write(6,*)' *** Warning: arguments to DFMD must be nonnegative!'
        write(6,*)' ***          I will take their ABS !!! '
        dFMD=abs(dFMD)
        Xsigma=abs(Xsigma)
       EndIf
      EndIf
*---  Process BKAP command for BK type of approximation (Giovanni Li Manni J.:GLMJ) Nov 2011-------------*
*
      If (KeyBKAP) Then
       DoBKAP = .true.
       Call SetPos_m(LUInput,'BKAP',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after BKAP keyword.'
       Read(LUInput,*,End=9910,Err=9920) NGASBK
       Read(LUInput,*,End=9910,Err=9920)(IOCCPSPC(IGAS,1),IGAS=1,NGASBK)
       Read(LUInput,*,End=9910,Err=9920)(IOCCPSPC(IGAS,2),IGAS=1,NGASBK)
       ReadStatus=' O.K. reading data after BKAP keyword.'
       If (DBG) Then
         write(6,*) ' BKAP: BK-type of approximation in action'
        write(6,*) ' Min and Max for subspace with exact Hamiltonian '
        write(6,*) ' =============================================== '
        write(6,*) ' NGASBK :',NGASBK
        write(6,*) '              Min. Occ.      Max. Occ.           '
         Do IGAS = 1, NGASBK
           write(6,'(A,I2,10X,I3,9X,I3)')
     &     '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
         End do
       End If
      End If
*
*---  Process SPLI command for SplitCAS calculations (Giovanni Li Manni J.:GLMJ) -------------*
*
      If (KeySPLI) Then
       If (DBG) Then
         Write(6,*) ' SPLI (Activation SplitCAS)'
       End If
       DoSplitCAS = .true.
       EnerSplit  = .true.
       iDimBlockA = 0
**** The energy gap (GapSpli) is in eV ****
       GapSpli    = 13.61
       lrootSplit= 1
       thrSplit = 1.0d-6
       MxIterSplit = 50
       Call SetPos_m(LUInput,'SPLI',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
      End If
*
*---  Process NUSP command for Numerical SplitCAS param. (GLMJ) ------------*
*
      If (KeyNUSP) Then
       If (DBG) Then
       Write(6,*)' NUSP - Manual Setting of Numerical SplitCAS Param.'
       End If
       EnerSplit  = .false.
       NumSplit  = .true.
       Call SetPos_m(LUInput,'NUSP',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after NUSP keyword.'
       Read(LUInput,*,End=9910,Err=9920)
     *              lrootSplit,iDimBlockA,MxIterSplit
       ReadStatus=' O.K. reading data after NUSP keyword.'
*       Read(LUInput,*,End=9910,Err=9920) lrootSplit
*       Read(LUInput,*,End=9910,Err=9920) iDimBlockA
*       Read(LUInput,*,End=9910,Err=9920) MxIterSplit
       ReadStatus=' Failure reading data after NUSP keyword.'
       Read(LUInput,*,End=9910,Err=9920) ThrSplit
       ReadStatus=' O.K. reading data after NUSP keyword.'
       If (DBG) then
         Write(6,*) ' Root to be opt. in SplitCAS = ', lrootSplit
         Write(6,*) ' AA block size in SplitCAS = ',iDimBlockA
         Write(6,*) ' Max iteration in SplitCAS = ', MxIterSplit
         Write(6,*) ' Root to be opt. in SplitCAS = ', ThrSplit
       end if
      End If
*
*---  Process ENSP command for Energetical SplitCAS param. (GLMJ) ------------*
*
      If (KeyENSP) Then
       If (DBG) Then
         Write(6,*)
     *   ' ENSP - Manual Setting of Energetical SplitCAS Param.'
       End If
       EnerSplit  = .true.
       NumSplit  = .false.
       iDimBlockA = 0
       Call SetPos_m(LUInput,'ENSP',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after ENSP keyword.'
       Read(LUInput,*,End=9910,Err=9920)
     *              lrootSplit,GapSpli,MxIterSplit
       ReadStatus=' O.K. reading data after ENSP keyword.'
*       Read(LUInput,*,End=9910,Err=9920) lrootSplit
*       Read(LUInput,*,End=9910,Err=9920) iDimBlockA
*       Read(LUInput,*,End=9910,Err=9920) MxIterSplit
       ReadStatus=' Failure reading data after ENSP keyword.'
       Read(LUInput,*,End=9910,Err=9920) ThrSplit
       ReadStatus=' O.K. reading data after ENSP keyword.'
      If (DBG) Write(6,*)' Root to be opt. in SplitCAS = ', lrootSplit
       If (DBG) Write(6,*) ' Energy gap in SplitCAS = ', GapSpli
      If (DBG) Write(6,*) ' Max iteration in SplitCAS = ', MxIterSplit
       If (DBG) Write(6,*) ' Root to be opt. in SplitCAS = ', ThrSplit
      End If
*
*---  Process PESP command for Percentage SplitCAS param. (GLMJ) -------*
*
      If (KeyPESP) Then
       If (DBG) Then
         Write(6,*)
     *   ' PESP - Manual Setting of Percentage SplitCAS Param.'
       End If
       EnerSplit  = .false.
       PerSplit  = .true.
       NumSplit  = .false.
       iDimBlockA = 0
       Call SetPos_m(LUInput,'PESP',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after PESP keyword.'
       Read(LUInput,*,End=9910,Err=9920)
     *              lrootSplit,PercSpli,MxIterSplit
       Read(LUInput,*,End=9910,Err=9920) ThrSplit
       ReadStatus=' O.K. reading data after PESP keyword.'
      If (DBG) Write(6,*) ' Root to be opt. in SplitCAS = ',lrootSplit
       If (DBG) Write(6,*) ' Percentage in SplitCAS = ', PercSpli
      If (DBG) Write(6,*) ' Max iteration in SplitCAS = ', MxIterSplit
       If (DBG) Write(6,*) ' Root to be opt. in SplitCAS = ', ThrSplit
      End If
*
*------- Process FOSP command for First Order SplitCAS Approx. (GLMJ)  -------*
*
      If (KeyFOSP) Then
       If (DBG) Then
         Write(6,*)
     *   ' FOSP - First Order SplitCAS Approx.'
       End If
       FOrdSplit  = .true.
      End If
*
*---  Process OPTO keyword: Optimal Output for RASSCF/CASPT2
*                           optimizations - GG Nov 2008.
      If (KeyOPTO) Then
       If (DBG) Then
        Write(6,*) ' OPTO keyword was used.'
        Write(6,*) '(Optimal Output for RASSCF/CASPT2)'
       End If
       lOPTO=.True.
       Call SetPos_m(LUInput,'OPTO',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process SXDAmp command --------------------------------------------*
      If (KeySXDA) Then
       If (DBG) Write(6,*)' SXDAMPING was requested.'
       Call SetPos_m(LUInput,'SXDA',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after SXDAMP keyword.'
       Read(LUInput,*,End=9910,Err=9920) SXDamp
       ReadStatus=' O.K. after reading data after SXDAMP keyword.'
       If (DBG) Write(6,*)' Parameter SXDamp=',SXDamp
       Call ChkIfKey_m()
      End If

*---  Process LOWM command ---
      If (KeyLOWM) Then
       If (DBG) Then
         Write(6,*)' LOWM keyword was used to force the CI routines'
         Write(6,*)' to use Slater Determinants for low M and not M=S'
       End If
       LOWMS=1
       Call SetPos_m(LUInput,'LOWM',Line,iRc)
       Call ChkIfKey_m()
      End If

*---  Process LOWD keyword: Turn on Lowdin orthonormalization of CMOs
      If (KeyLOWD) Then
       If (DBG) Then
         Write(6,*)' LOWDIN orthonormalization was requested'
         Write(6,*)' but from Jan 12 2010 that is default anyway.'
       End If
       Lowdin_ON=.True.
       Call SetPos_m(LUInput,'LOWD',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*
*---  Process PRWF command --------------------------------------------*
      If (KeyPRWF) Then
       If (DBG) Write(6,*)' The PRWF keyword was used.'
       Call SetPos_m(LUInput,'PRWF',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after PRWF keyword.'
       Read(LUInput,*,End=9910,Err=9920) PRWTHR
       ReadStatus=' O.K. after reading data after PRWF keyword.'
      If (DBG) Write(6,*)' Print CI coefficients larger than PRWTHR=',
     &                     PRWTHR
       Call ChkIfKey_m()
      End If
*
*
*---  Process PRSD command --------------------------------------------*
      If (KeyPRSD) Then
       If (DBG) Write(6,*)' The PRSD keyword was used.'
       Call SetPos_m(LUInput,'PRSD',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       If (DBG) Write(6,*)' Print determinant expansions of CSFs'
       Call ChkIfKey_m()
      End If
*
*
*---  Process ALPH command --------------------------------------------*
      If (KeyALPH) Then
       If (DBG) Write(6,*)' The ALPH keyword was used.'
       Call SetPos_m(LUInput,'ALPH',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after ALPH keyword.'
       Read(LUInput,*,End=9910,Err=9920) iAlphaBeta
       ReadStatus=' O.K. after reading data after ALPH keyword.'
       If (iAlphaBeta.lt.0) iAlphaBeta=-1
       If (iAlphaBeta.gt.0) iAlphaBeta=1
       If (DBG) Then
        If (iAlphaBeta.eq.1) Write(6,*)' Read alpha orbitals from UHF'
        If (iAlphaBeta.eq.-1) Write(6,*)' Read beta orbitals from UHF'
       End If
       Call ChkIfKey_m()
      End If
*
*---  Process ALPH command --------------------------------------------*
      If (KeyFARO) Then
        DoFaro = .TRUE.
      End If
*
*---  Process NOCA command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if NOCALC case.'
      If (KeyNOCA) Then
       If (DBG) Write(6,*) ' NOCALC keyword was used.'
       INOCALC=1
       Call SetPos_m(LUInput,'NOCA',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process SAVE command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if SAVE_EXP case.'
      If (KeySAVE) Then
       If (DBG) Write(6,*) ' SAVE_EXP keyword was used.'
       ISAVE_EXP=1
       Call SetPos_m(LUInput,'SAVE',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process EXPA command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if EXPAND case.'
      If (KeyEXPA) Then
       If (DBG) Write(6,*) ' EXPAND keyword was used.'
       IEXPAND=1
       Call SetPos_m(LUInput,'EXPA',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process GRAD command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if GRADient case.'
      If (KeyGRAD) Then
       If (DBG) Write(6,*) ' GRADient keyword was used.'
       DoGradPDFT=.true.
       Call SetPos_m(LUInput,'GRAD',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process NOGR command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if NOGRadient case.'
      If (KeyNOGR) Then
       If (DBG) Write(6,*) ' NOGRadient keyword was used.'
       DoNoGrad=.true.
       Call SetPos_m(LUInput,'GRAD',Line,iRc)
       Call ChkIfKey_m()
      End If
*
*---  Process DMRG command --------------------------------------------*
#ifdef _ENABLE_BLOCK_DMRG_
      If (KeyDMRG) Then
* NN.14 FIXME: When DMRG option is disabled at compilation,
*       this should give an error, but just ignored for the time.
       If (DBG) Then
         Write(6,*) ' DMRG (Use DMRG algorithm instead of FCI)'
       End If
       Call SetPos_m(LUInput,'DMRG',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure reading data after DMRG keyword.'
       Read(LUInput,*,End=9910,Err=9920) MxDMRG
       ReadStatus=' O.K. after reading data after DMRG keyword.'
       If (DBG) Write(6,*) ' Nr. of states=',MxDMRG
       DoDMRG=.True.
       Call ChkIfKey_m()
      End If
*---  Process 3RDM command --------------------------------------------*
      If (Key3RDM) Then
       If (DBG) Then
         Write(6,*) ' 3RDM (Compute 3RDM for DMRG-Cu4-CASPT2)'
       End If
       Do3RDM=.True.
       Call SetPos_m(LUInput,'3RDM',Line,iRc)
       Call ChkIfKey_m()
      End If
#endif
*---  All keywords have been processed ------------------------------*

************************************************************************
* Generate artificial splitting or RAS into GAS for parallel blocking  *
************************************************************************
      IF (.NOT.IDOGAS) THEN
* SVC: convert CAS/RAS to general GAS description here, then we only
* need to copy it for lucia later, which always uses GAS description.
        NGSSH(1,1:NSYM)=NRS1(1:NSYM)
        NGSSH(2,1:NSYM)=NRS2(1:NSYM)
        NGSSH(3,1:NSYM)=NRS3(1:NSYM)
        IGSOCCX(1,1) = MAX(2*SUM(NRS1(1:NSYM))-NHOLE1,0)
        IGSOCCX(1,2) = 2*SUM(NRS1(1:NSYM))
        IGSOCCX(2,1) = NACTEL - NELEC3
        IGSOCCX(2,2) = NACTEL
        IGSOCCX(3,1) = NACTEL
        IGSOCCX(3,2) = NACTEL
      END IF
*
!Consideraations for gradients/geometry optimizations

*     Numerical gradients requested in GATEWAY
      Call Qpg_iScalar('DNG',DNG)
      If (DNG) Then
         Call Get_iScalar('DNG',iDNG)
         DNG = iDNG.eq.1
      End If
      DNG=DoNoGrad.or.DNG
*
*     Inside LAST_ENERGY we do not need analytical gradients
      ProgName=Get_SuperName()
      If (ProgName(1:11).eq.'last_energy') DNG=.true.
*
*     Inside NUMERICAL_GRADIENT override input!
      If (ProgName(1:18).eq.'numerical_gradient') Then
         DNG=.true.
         DoGradPDFT=.false.
      End If
*
*     Check to see if we are in a Do While loop
         Call GetEnvF('EMIL_InLoop',emiloop)
         If (emiloop.eq.' ') emiloop='0'
         Call GetEnvF('MOLCAS_IN_GEO',inGeo)
         If ((emiloop(1:1).ne.'0') .and. inGeo(1:1) .ne. 'Y'
     &       .and. .not.DNG) Then
            DoGradPDFT=.true.
         End If
*                                                                      *
************************************************************************
*                                                                      *
*     Select default root for geometry optimization
*
      If (NROOTS.gt.1.and.irlxroot.eq.0)  Then
*
*        Check if multi state SA-CASSCF
*
         nW = 0
         Do iR = 1, LROOTS
            If (WEIGHT(iR).ne.0.0D0) nW = nW + 1
         End Do
         If (nW.ne.1) Then
            iRlxRoot=iroot(LROOTS)
         Else
            Do iR = 1, LROOTS
               If (WEIGHT(iR).ne.0.0D0) iRlxRoot=iroot(iR)
            End Do
         End If
      End If
      If (NROOTS.eq.1.or.LROOTS.eq.1) iRlxRoot=iRoot(1)
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute IZROT. IZROT is a matrix (lower triangular over the -----*
*     active space), which specifies which t,u rotations should be
*     avoided, since the orbitals belong to the same RAS space.
*     This is the only way the RAS concept is explicitly used in the
*     SX section of the program.
      ITU=0
      DO ISYM=1,NSYM
        NAO=NASH(ISYM)
*
        IF(DOBKAP)THEN
*.Giovanni... BK stuff SplitCAS related. We want to treat RAS CI space as CAS.
          DO NT=2,NAO
            DO NU=1,NT-1
              ITU=ITU+1
              IZROT(ITU)=1
            END DO
          END DO
        ELSE
*
*
*
          DO NT=2,NAO
            DO NU=1,NT-1
              ITU=ITU+1
              IZROT(ITU)=0
CSVC: check if NU<NT are included in the same gas space
              NGSSH_LO=0
              DO IGAS=1,NGAS
                NGSSH_HI=NGSSH_LO+NGSSH(IGAS,ISYM)
                IF (NU.GT.NGSSH_LO.AND.NT.LE.NGSSH_HI) THEN
                  IZROT(ITU)=1
                END IF
                NGSSH_LO=NGSSH_HI
              END DO
            END DO
          END DO
        END IF
      END DO
*
      Call Put_iArray('nIsh',nIsh,nSym)
      Call Put_iArray('nAsh',nAsh,nSym)
      Call Put_iScalar('Multiplicity',ISPIN)
*
*---  Initialize Cholesky information if requested
      if (DoCholesky)then
         Call Cho_X_init(irc,ChFracMem)
         if (irc.ne.0) Go To 9930
      endif

* ===============================================================

 100  CONTINUE
*PAM Jump here in case of CASVB (IFVB.eq.2)
      If (DBG) Write(6,*)' After 100 CONTINUE.'

*PAM July 2007 Check in case of CI restart:
      If (DBG) Write(6,*)' Check if CI-Restart.'
      If (KeyCIRE) Then
* Test read:
        If (DBG) Write(6,*)' Yes it is!'
        iJOB=-1
        Call f_Inquire('JOBOLD',lExists)
        If (lExists) Then
           If (DBG) Write(6,*)' ''JOBOLD'' exists.'
           iJOB=1
        Else
           Call f_Inquire(IPHNAME,lExists)
           If (lExists) Then
             iJOB=0
             If (DBG) Write(6,*)' No ''JOBOLD'', but JOBIPH exists.'
             If (DBG) Write(6,*)' It is named ',IPHNAME
           End If
        End If
        If (iJOB.eq.1.or.iJOB.eq.0) Then
           if(JOBOLD.le.0) Then
             JOBOLD=20
           end if
           If (iJOB.eq.1) Then
              Call DaName(JOBOLD,'JOBOLD')
           Else
              Call DaName(JOBOLD,IPHNAME)
           End If
           IDisk=0
           Call IDaFile(JOBOLD,99,IScratch,10,IDisk)
           If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
             Call DaClos(JOBOLD)
             JOBOLD=-1
           Else If (JOBOLD.gt.0) Then
             JOBOLD=-1
           End If
        Else
           IScratch(1)=0
        End If
        If (IScratch(1).eq.0) Then
C Test read failed. JOBOLD cannot be used.
        Write(6,*)' Test read shows that there is no usable interface'
          Write(6,*)' file, necessary for the requested CI restart.'
          Write(6,*)' Most probable reason: the user has forgotten to'
          Write(6,*)' provide this file. The program will continue,'
          Write(6,*)' but there can be no CI restart.'
          ICIRST=0
        End If
      End If
*PAM July 2007 End of addition
* ===============================================================
*
*     Initialize seward
*
      If (DBG) Write(6,*)' Initialize seward.'
      nDiff = 0
      If (DSCF           .or.
     &    RF_On()        .or.
     &    Langevin_On()  .or.
     &    PCM_On()       .or.
     &    Do_OFEmb       .or.
     &    KSDFT.ne.'SCF'     )
     &    Call IniSew(Info,DSCF.or.Langevin_On().or.PCM_On(),nDiff)
* ===============================================================
*
*     Check the input data
*
      If (DBG) Then
        Write(6,*)' Call ChkInp.'
        Call XFlush(6)
      End If
      Call ChkInp_m()
* ===============================================================
*
*     In DMRG-CASSCF, skip GUGA and LUCIA settings
*
      NCONF=1
      GoTo 9000
      If(DoDMRG) GoTo 9000
* ===============================================================
*
*     Construct the Guga tables
*
      if(.not.iDoGas) then
        call gugactl_m
      else
        call mknsm_m
      end if
* ===============================================================
*
*     Construct the determinant tables
*

      If (DBG) Write(6,*)' Construct the determinant tables.'
      MS2 = iSpin - 1
*
* Set variables needed in Lucia_Ini
*
      Call iCopy(mxGAS*mxSym,ngssh,1,ngssh_Molcas,1)
      Call iCopy(mxGAS*2,igsoccx,1,igsoccx_Molcas,1)
      Call iCopy(nSym,norb,1,norb_Molcas,1)
      Call iCopy(nSym,nbas,1,nbas_Molcas,1)
      Call iCopy(nSym,nish,1,nish_Molcas,1)
      potnuc_Molcas    = potnuc
      thre_Molcas      = thre
      nsym_Molcas      = nsym
      nactel_Molcas    = nactel
      ms2_Molcas       = ms2
      ispin_Molcas     = ispin
      lsym_Molcas      = lsym
      NHOLE1_Molcas    = NHOLE1
      NELEC3_Molcas    = NELEC3
      itmax_Molcas     = itmax
      rtoi_Molcas      = rtoi
      nroots_Molcas    = Max(nroots,lRoots)
      ipt2_Molcas      = ipt2
      iprci_molcas     = iprloc(3)
      ngas_molcas      = ngas
      INOCALC_MOLCAS   = INOCALC
      ISAVE_EXP_MOLCAS = ISAVE_EXP
      IEXPAND_MOLCAS   = IEXPAND

*
* And call Lucia_Ini to initialize LUCIA
*
* Combinations don't work for CASVB (at least yet)!
      If (ifvb .ne. 0) iSpeed(1) = 0
*
      CALL Lucia_Util('Ini',iDummy,iDummy,Dummy)
* to get number of CSFs for GAS
      nconf=0
      do i=1,mxsym
        nconf=nconf+ncsasm(i)
      end do
*
      ISCF=0
      IF (ISPIN.EQ.NAC+1.AND.NACTEL.EQ.NAC) ISCF=1
      IF (ISPIN.EQ.1.AND.NACTEL.EQ.2*NAC)   ISCF=1
      IF (ISCF.EQ.1) THEN
         NCONF=1
         MAXJT=1
      END IF
*
*     If the CI-root selectioning option has been specified translate
*     the reference configuration numbers from the split graph GUGA
*     to the symmetric group numbering
*
* ===============================================================
      IF (ICICH.EQ.1) THEN
        CALL GETMEM('UG2SG','ALLO','INTE',LUG2SG,NCONF)
        CALL UG2SG_m(NROOTS,NCONF,NAC,NACTEL,LSYM,IPR,
     *             IWORK(KICONF(1)),IWORK(KCFTP),IWORK(LUG2SG),
     *             ICI,JCJ,CCI,MXROOT)
        CALL GETMEM('UG2SG','FREE','INTE',LUG2SG,NCONF)
      END IF
* ===============================================================

#ifdef module_DMRG
!     call set_dmrg_cfg()
#endif

      ! faroald initializations
!      DOFARO = .True.
      IF (DOFARO) THEN
        IF (NSYM.GT.1) THEN
          WRITE(6,'(1X,A)') 'FARO keyword was used, but NSYM > 1,'
          WRITE(6,'(1X,A)') 'switching to LUCIA as the CI backend.'
          DOFARO=.FALSE.
        ELSE
          WRITE(6,'(1X,A)') '**EXPERIMENTAL**'
          WRITE(6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
          WRITE(6,'(1X,A)') '**EXPERIMENTAL**'
          CALL FAROALD_INIT(NACTEL,NASH(1),ISPIN)
          CALL CITRANS_INIT(NACTEL,NASH(1),ISPIN)
        END IF
      END IF
      DOFARO = .FALSE.

      Go to 9000
*
*---  Error exits -----------------------------------------------------*
9810  CONTINUE
      If (IPRLEV.ge.TERSE) Then
       Call WarningMessage(2,'Error in input preprocessing.')
       Write(6,*)' PROC_INP: A keyword was found during prescanning'
       Write(6,*)' the input file, but when later trying to locate'
       Write(6,*)' this input, it could not be found. Something has'
       Write(6,*)' happened to the input file, or else there is some'
       Write(6,*)' strange program error.'
       iRc=_RC_INPUT_ERROR_
      End If
      Go to 9900

9910  CONTINUE
      Call WarningMessage(2,'End of input file during preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9910 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900
*
9920  CONTINUE
      Call WarningMessage(2,'Read error during input preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9920 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900
*
9930  CONTINUE
      Call WarningMessage(2,'Error during input preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9930 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900

*---  Normal exit -----------------------------------------------------*
9000  CONTINUE
      close(989)
      If (DBG) Write(6,*)' Normal exit from PROC_INP.'
      Call qExit('Proc_Inp')
      Return
*---  Abnormal exit -----------------------------------------------------*
9900  CONTINUE
      If (DBG) Write(6,*)' Abnormal exit from PROC_INP.'
      Call qExit('Proc_Inp')
      Return
      End
