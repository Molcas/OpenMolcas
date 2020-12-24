************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Markus P. Fuelscher                              *
*               1995, Martin Schuetz                                   *
*               2004,2005, Thomas Bondo Pedersen                       *
************************************************************************
      Subroutine RdInp(CMO,Eall,Eocc,Eext,iTst,ESCF)
************************************************************************
*                                                                      *
*     Locate input stream and read commands                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     modified by MGS                                                  *
*                 TBP                                                  *
*     University of Lund, Sweden, 1992/95, 2004/05.                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Dimension CMO(*),Eall(*),Eocc(*),Eext(*)
#include "cddos.fh"
#include "chomp2_cfg.fh"
#include "real.fh"
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "WrkSpc.fh"
#include "print_mbpt2.fh"
#include "warnings.fh"
#include "Molcas.fh"
#include "namact.fh"
      Parameter ( nCom=43 )
      Character*4 Command,ComTab(nCom)
      Character*180 Get_Ln
      External Get_Ln
      Character*100 ProgName, Get_SuperName
      External Get_SuperName
      Integer  iPrintLevel
      External iPrintLevel
      Logical  ChoMP2_ChkPar
      External ChoMP2_ChkPar
      Logical  Reduce_Prt
      External Reduce_Prt
      Data ComTab/'TITL','FROZ','DELE','SFRO','SDEL',
     &            'EXTR','PRIN','TEST','TSTP','PRPT',
     &            'LUMO','EREF','VIRA','T1AM','GRDT',
     &            'LAPL','GRID','BLOC','CHOA','$$$$',
     &            '$$$$','$$$$','DECO','NODE','THRC',
     &            'SPAN','MXQU','PRES','CHKI','FORC',
     &            'VERB','NOVE','FREE','PREC','SOSM',
     &            'OEDT','OSFA','LOVM','DOMP','FNOM',
     &            'GHOS','NOGR','END '/
      Character*180 Line
      Character*8 emiloop
      Character*8 inGeo
      Character*72  Blank
      Logical FrePrt,ERef_UsrDef, DecoMP2_UsrDef,DNG,NoGrdt
      Logical lTit,lFro,lFre,lDel,lSFro,lSDel,lExt,lPrt,LumOrb
      Character*80 VecTitle
      Real*8 ESCF
      Dimension iDummy(1)
*----------------------------------------------------------------------*
*     Locate "start of input"                                          *
*----------------------------------------------------------------------*
      lTit=.false.
      lFro=.false.
      lFre=.false.
      lDel=.false.
      lSFro=.false.
      lSDel=.false.
      lExt=.false.
      lPrt=.false.
      LovMP2=.false.
      DoMP2=.false.
      FNOMP2=.false.
      LumOrb=.false.
      all_vir=.false.
      DoT1amp=.false.
      ESCF=0.0d0
      Thr_ghs=5.0d-1
      DelGHost=.false.
      ERef_UsrDef=.false.
      DecoMP2_UsrDef=.false.
      Call Put_iScalar('mp2prpt',0)
*
*     copy input from standard input to a local scratch file
*
      LuSpool = 17
      Call SpoolInp(LuSpool)
      LuRd=5
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'MBPT2')
      Blank=' '
*----------------------------------------------------------------------*
*     Define default values                                            *
*----------------------------------------------------------------------*
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0



*
      If (DoCholesky) Then
         ChoAlg=2
      Else
         ChoAlg=-999999
      End If
      DecoMP2=Decom_Def
      ThrMP2=-9.9D9
      SpanMP2=Span_Def
      MxQualMP2=MxQual_Def
      ChkDecoMP2=.false.
      ForceBatch=.false.
      If (iPL .ge. 3) Then
         Verbose=.true.
      Else
         Verbose=.false.
      End If
* --- Scaled Opposite-Spin MP2
      SOS_mp2=.false.
      set_cd_thr=.true.
      OED_Thr=1.0d-8
      C_os=1.3d0
      EOSMP2=0.0d0
* --- Frozen natural orbitals
      DoFNO=.false.
* --- MP2-gradient/1pdm
      DoDens=.false.
      DoGrdt=.false.
      NoGrdt=.false.
      NoGamma=.false.
* --- Laplace transform
      Laplace=.false.
      Laplace_nGridPoints=0
      Laplace_BlockSize=Laplace_BlockSize_Def
* --- LDF settings
      If (DoLDF) Then
         SOS_MP2=.true.
         Laplace=.true.
      End If

* ---------------------
      nTit=0
      iTst=0
C     nTstP=-1
      iPrt=0
      Call Get_iArray('Non valence orbitals',nFro1,nSym)
      Do iSym=1,nSym
C        nFro1(iSym)=0
         nFro2(iSym)=0
         Do iOrb=1,mxFro
            iFro(iSym,iOrb)=0
         End Do
         nDel1(iSym)=0
         nDel2(iSym)=0
         Do iOrb=1,mxFro
            iDel(iSym,iOrb)=0
         End Do
      End Do
      nFre=0
*----------------------------------------------------------------------*
*     Read the input stream line by line and identify key command      *
*----------------------------------------------------------------------*
100   Line = Get_Ln(LuSpool)
      Call StdFmt(Line,Command)
      jCom=0
      Do iCom=1,nCom
        If ( Command.eq.ComTab(iCom) ) jCom=iCom
      End Do
      If ( jCom.eq.0 ) Then
         Write (6,*) 'RdInp: Illegal keyword!'
         Write (6,'(A,A)') 'Command=',Command
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Branch to the processing of the command sections                 *
*----------------------------------------------------------------------*
 110  Continue
      Select Case(jCom)
         Case(1)
            Goto 501
         Case(2)
            Goto 502
         Case(3)
            Goto 503
         Case(4)
            Goto 504
         Case(5)
            Goto 505
         Case(6)
            Goto 506
         Case(7)
            Goto 507
         Case(8)
            Goto 508
C        Case(9)
C           Goto 509
         Case(10)
            Goto 510
         Case(11)
            Goto 511
         Case(12)
            Goto 512
         Case(13)
            Goto 513
         Case(14)
            Goto 514
         Case(15)
            Goto 515
         Case(16)
            Goto 516
         Case(17)
            Goto 517
         Case(18)
            Goto 518
         Case(19)
            Goto 519
         Case(20)
            Goto 520
         Case(21)
            Goto 521
         Case(22)
            Goto 522
         Case(23)
            Goto 523
         Case(24)
            Goto 524
         Case(25)
            Goto 525
         Case(26)
            Goto 526
         Case(27)
            Goto 527
         Case(28)
            Goto 528
         Case(29)
            Goto 529
         Case(30)
            Goto 530
         Case(31)
            Goto 531
         Case(32)
            Goto 532
         Case(33)
            Goto 533
         Case(34)
            Goto 534
         Case(35)
            Goto 535
         Case(36)
            Goto 536
         Case(37)
            Goto 537
         Case(38)
            Goto 538
         Case(39)
            Goto 539
         Case(40)
            Goto 540
         Case(41)
            Goto 541
         Case(42)
            Goto 542
         Case Default
            Goto 1000
      End Select
*---  Process the "TITL" input card -----------------------------------*
 501  Continue
      If ( lTit ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'Title option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lTit=.true.
 15   Line = Get_Ln(LuSpool)
      Call StdFmt(Line,Command)
      jCom=0
      Do iCom=1,nCom
        If ( Command.eq.ComTab(iCom) ) jCom=iCom
      End Do
      If ( jCom.ne.0 ) Goto 110
      nTit=nTit+1
      If ( nTit.le.mxTit ) then
         Title(nTit)=Line(1:80)
         Goto 15
      End If
      Goto 100
*---  Process the "FROZ" input card -----------------------------------*
 502  Continue
      If ( lFre .or. lFro) Then
         Write (6,*) 'RdInp: Error while reading input!'
         If (lFro) Write (6,*) 'Frozen option already processed!'
         If (lFre) Write (6,*) 'Freeze option and Frozen option ',
     &                         'are incompatible!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lFro=.true.
      Line = Get_Ln(LuSpool)
*
      Write (6,*)
      Write (6,'(A)') 'WARNING!'
      Write (6,'(A)') 'Default frozen orbitals as non valence orbitals'
     &              //' is overwritten by user input.'
      Write (6,'(A,8I4)') 'Default values:',(nFro1(iSym),iSym=1,nSym)
      Write (6,*)
*
      Read(Line,*,err=995) (nFro1(iSym),iSym=1,nSym)
      Goto 100
*---  Process the "DELE" input card -----------------------------------*
 503  Continue
      If ( lDel ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'Delete option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lDel=.true.
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) (nDel1(iSym),iSym=1,nSym)
      Goto 100
*---  Process the "SFRO" input card -----------------------------------*
 504  Continue
      If ( lSFro ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'SFrozen option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lSFro=.true.
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) (nFro2(iSym),iSym=1,nSym)
      Do iSym=1,nSym
         Line = Get_Ln(LuSpool)
         Read(Line,*,err=995)
     &       (iFro(iSym,iOrb),iOrb=1,nFro2(iSym))
      End Do
      Goto 100
*---  Process the "SDEL" input card -----------------------------------*
 505  Continue
      If ( lSDel ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'SDelete option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lSDel=.true.
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) (nDel2(iSym),iSym=1,nSym)
      Do iSym=1,nSym
         Line = Get_Ln(LuSpool)
         Read(Line,*,err=995)
     &       (iDel(iSym,iOrb),iOrb=1,nDel2(iSym))
      End Do
      Goto 100
*---  Process the "Extract" input card --------------------------------*
 506  Continue
      If ( lExt ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'Extract option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      Write (6,*) 'RdInp: EXTRACT option is redundant and is ignored!'
      Goto 100
*---  Process the "Print" input card ----------------------------------*
 507  Continue
      If ( lPrt ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'Print option already processed!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
      lPrt=.true.
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) iPrt
      Goto 100
*---  Process the "Test" input card -----------------------------------*
 508  Continue
      iTst=1
      Goto 100
*---  Process the "TSTP" input card -----------------------------------*
C509  Continue
C     iTst=1
C     Line = Get_Ln(LuSpool)
C     Read(Line,*,err=995) nTstP
      Goto 100
*---  Process the "PRPT" input card -----------------------------------*
 510  Continue
      Call Put_iScalar('mp2prpt',1)
      DoDens = .true.
      NoGamma = .true.
      Goto 100
*---  Process the "LUMO" input card -----------------------------------*
 511  Continue
      LumOrb = .true.
      Goto 100
*---  Process the "EREF" input card -----------------------------------*
 512  Continue
      If (.not.LumOrb) Then
         Write (6,*) 'RdInp: Error while reading input!'
         Write (6,*) 'EREF can be used only with LumOrb.'
         Write (6,*) '(Note: LumOrb keyword must precede EREF)'
         Call Abend()
      EndIf
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) ESCF
      ERef_UsrDef=.true.
      Goto 100
*---  Process the "VIRA" input card -----------------------------------*
 513  Continue
      all_Vir=.true.
      Goto 100
*---  Process the "T1AM" input card -----------------------------------*
 514  Continue
      DoT1amp=.true.
      If (.not.DoCholesky) Then
         Write (6,*) 'RdInp: T1AM is available only with Cholesky/RI .'
         Call Abend()
      EndIf
      Goto 100
*---  Process the "GRDT" input card -----------------------------------*
 515  Continue
      Call Put_iScalar('mp2prpt',2)
      DoDens = .true.
      DoGrdt = .true.
      Goto 100
*---  Process the "LAPLace" input card --------------------------------*
*     Laplace transform with default or previously specified grid.
 516  Continue
      Laplace=.true.
      Goto 100
*---  Process the "GRID" input card -----------------------------------*
*     Read number of Laplace grid points (activates Laplace as well)
 517  Continue
      Laplace=.true.
      Line=Get_Ln(LuSpool)
      Read(Line,*,Err=995) Laplace_nGridPoints
      Laplace_nGridPoints=max(0,Laplace_nGridPoints)
      If (Laplace_nGridPoints.gt.Laplace_mGridPoints) Then
         Call WarningMessage(2,'Input Error')
         Write(6,'(A,I6)') 'Number of Laplace grid points specified:',
     &                     Laplace_nGridPoints
         Write(6,'(A,I6)') 'Max allowed:                            ',
     &                     Laplace_mGridPoints
         Call Quit(_RC_INPUT_ERROR_)
      End If
      Goto 100
*---  Process the "BLOC" input card -----------------------------------*
*     Read vector block size for CD/DF-Laplace-SOS-MP2
*     (Activates Laplace as well)
 518  Continue
      Laplace=.true.
      Line=Get_Ln(LuSpool)
      Read(Line,*,Err=995) Laplace_BlockSize
      Laplace_BlockSize=max(0,Laplace_BlockSize)
      If (Laplace_BlockSize.eq.0) Then
         Laplace_BlockSize=Laplace_BlockSize_Def
      End If
      GoTo 100
*---  Process the "CHOAlgorithm" input card ---------------------------*
 519  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,Err=995) ChoAlg
      If (ChoAlg .lt. 0) Then
         ChoAlg = 0
      Else If (ChoAlg .gt. 2) Then
         ChoAlg = 2
      End If
      Goto 100
*---  Process the "$$$$" input card -----------------------------------*
 520  Continue
c     not used
      Goto 100
*---  Process the "$$$$" input card -----------------------------------*
 521  Continue
c     not used
      Goto 100
*---  Process the "$$$$" input card -----------------------------------*
 522  Continue
c     not used
      Goto 100
*---  Process the "DECOmpose MP2 integrals" card ----------------------*
 523  Continue
      DecoMP2=.true.
      DecoMP2_UsrDef=.true.
      Goto 100
*---  Process the "NODEcomposition of MP2 integrals" card -------------*
 524  Continue
      DecoMP2=.false.
      DecoMP2_UsrDef=.true.
      Goto 100
*---  Process the "THRCholesky" card ----------------------------------*
 525  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) ThrMP2
      set_cd_thr=.false.
      Goto 100
*---  Process the "SPAN" card -----------------------------------------*
 526  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) SpanMP2
      Goto 100
*---  Process the "MXQUal" input card ---------------------------------*
 527  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) MxQualMP2
      If (MxQualMP2 .lt. 1) MxQualMP2 = MxQual_Def
      Goto 100
*---  Process the "PRESort input card (OBSOLETE) ----------------------*
 528  Continue
      Goto 100
*---  Process the "CHKI" card -----------------------------------------*
 529  Continue
      ChkDecoMP2=.true.
      Goto 100
*---  Process the "FORCebatching" card --------------------------------*
 530  Continue
      ForceBatch=.true.
      Goto 100
*---  Process the "VERBose" card --------------------------------------*
 531  Continue
      Verbose=.true.
      Goto 100
*---  Process the "NOVErbose" card ------------------------------------*
 532  Continue
      Verbose=.false.
      Goto 100
*---  Process the "FREEze" card  --------------------------------------*
 533  Continue
      If ( lFre .or. lFro ) Then
         Write (6,*) 'RdInp: Error while reading input!'
         If (lFre) Write (6,*) 'Freeze option already processed!'
         If (lFro) Write (6,*) 'Frozen option and Freeze option ',
     &                         'are incompatible!'
         Write (6,'(A,A)') 'Last read line:',Line
         Call Abend()
      End If
*
      Write (6,*)
      Write (6,'(A)') 'WARNING!'
      Write (6,'(A)') 'Default frozen orbitals as non valence orbitals'
     &              //' is overwritten by user input.'
      Write (6,'(A,8I4)') 'Default values:',(nFro1(iSym),iSym=1,nSym)
      Write (6,*)
      Call ICopy(nSym,[0],0,nFro1,1)
*
      lFre=.true.
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) nFre
      Goto 100
*---  Process the "PRECision" card (= "THRCholesky" card) -------------*
 534  Continue
      Goto 525
* Commented out this line because it is not really possible to get here
* //Jonas Bostrom
*      Goto 100
*---  Process the "SOSMp2" card ---------------------------------------*
 535  Continue
      SOS_MP2=.true.
      If (.not.DoLDF) Then
         DecoMP2=.true.
         If (ChoMP2_ChkPar()) Then
            Call WarningMessage(2,'SOS-MP2 is not implemented '
     &                               //'for parallel runs. !! SORRY !!')
            Call Quit(_RC_NOT_AVAILABLE_)
         End If
      End If
      Goto 100
*---  Process the "OEDThreshold" card ---------------------------------*
 536  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) OED_Thr
      Goto 100
*---  Process the "OSFActor" card -------------------------------------*
 537  Continue
      Line = Get_Ln(LuSpool)
      Read(Line,*,err=995) C_os
      Goto 100
*---  Process the LovMP2 input ----------------------------------------*
 538  Continue
      If (.not. DoCholesky) Then
         WRITE(6,*)
         WRITE(6,*)'********************* ERROR ***********************'
         WRITE(6,*)' LovMP2 not implemented with conventional ERIs.'
         WRITE(6,*)' Please, use Cholesky or RI options.'
         WRITE(6,*)'***************************************************'
         Call Abend()
      Else
         LovMP2=.true.
      EndIf
      Read(LuSpool,*) nActa,ThrLov
*     nActa = number of active atoms
*     ThrLov = threshold for orbital selection
*
      If (ThrLov.lt.0.0d0 .or. ThrLov.ge.1.0d0) Then
         write(6,*)' Threshold out of range! Must be in [0,1[ '
         Call Abend()
      EndIf
*     namAct = names of active atoms (symm. indep. centers)
1538  Read(LuSpool,'(A)',End=995) Line
      If ( Line(1:1).eq.'*' ) Goto 1538
      If ( Line.eq.Blank) Goto 1538
      Call UpCase(Line)
      Do i=1,nActa
       Call LeftAd(Line)
       If( Line.eq.Blank) Goto 995
       j=Index(Line,' ')
       namAct(i)=Line(1:j-1)
       Line(1:j-1)=Blank(1:j-1)
      Enddo
      Goto 100
*
*---  Process DoMP2 card ----------------------------------------------*
 539  Continue
      DoMP2=.true.
      Goto 100
*
*---  process FNOM command
 540  Continue
      If (.not.DoCholesky) Then
         WRITE(6,*)
         WRITE(6,*)'********************* ERROR ***********************'
         WRITE(6,*)' FNO-MP2 not implemented with conventional ERIs.   '
         WRITE(6,*)' Please, use Cholesky or RI options.'
         WRITE(6,*)'***************************************************'
         Call Abend()
      Else
         FnoMP2=.true.
      EndIf
      Read(LuSpool,*) vkept
*
      If (vkept.le.0.0d0 .or. vkept.gt.1.0d0) Then
         write(6,*)' Requested fraction of virtual space out of range! '
         write(6,*)' Must be in ]0,1] '
         Call Abend()
      EndIf
      Goto 100
*
*>>>>>>>>>>>>> GHOST <<<<<< Removal of GHOST virtual space <<<<<<<<<
*
 541  Continue
      Read(LuSpool,'(A)',End=995) Line
      If ( Line(1:1).eq.'*' ) Goto 541
      If ( Line.eq.Blank ) Goto 541
      Read(Line,*,Err=995) Thr_ghs
      If (thr_ghs.lt.0.0d0 .or. thr_ghs.ge.1.0d0) Then
         write(6,*)' GHOST threshold out of range! Must be in [0,1[ '
         Call Abend()
      EndIf
      DelGHost=.true.
      Goto 100
*
*---  Process the "NOGR" input card -----------------------------------*
 542  Continue
      NoGrdt = .true.
      Goto 100
*
*---  Unused ----------------------------------------------------------*
c543  Continue
c     Goto 100
*----------------------------------------------------------------------*
*     "End of input"                                                   *
*----------------------------------------------------------------------*
1000  Continue
*
*---  Postprocessing for SOS-MP2 and Laplace
      If (SOS_MP2) Then
         If (.not.(DoCholesky.or.DoDF.or.DoLDF)) Then
            Call WarningMessage(2,
     &                         'SOS-MP2 only implemented for CD/DF/LDF')
            Call Quit(_RC_INPUT_ERROR_)
         End If
      End If
      If (Laplace) Then
         If (.not.(DoCholesky.or.DoDF.or.DoLDF)) Then
            Call WarningMessage(2,
     &          'Laplace transformation only implemented for CD/DF/LDF')
            Call Quit(_RC_INPUT_ERROR_)
         End If
         If (.not.SOS_MP2) Then
            Call WarningMessage(2,
     &  'Laplace transformation only implemented for CD/DF/LDF-SOS-MP2')
            Call Quit(_RC_INPUT_ERROR_)
         End If
         If (.not.DecoMP2_UsrDef) Then
            ! SOS-MP2 without Laplace automatically turns on DecoMP2.
            ! With Laplace, however, we should use the same default
            ! as for standard MP2. So, change back to default unless
            ! the user explicitly asked for something else.
            DecoMP2=Decom_Def
         End If
         If (DoDens .or. DoGrdt) Then
            Call WarningMessage(2,
     &          'Laplace transformation is incompatible with PRPT/GRDT')
            Call Quit(_RC_INPUT_ERROR_)
         End If
      End If
*
*     Check if the calculation is inside a loop and make analytical
*     gradients default in this case
*
*     Numerical gradients requested in GATEWAY
      Call Qpg_iScalar('DNG',DNG)
      If (DNG) Then
         Call Get_iScalar('DNG',iDNG)
         DNG = iDNG.eq.1
      End If
      DNG=NoGrdt.or.DNG
*
*     Inside LAST_ENERGY we do not need analytical gradients
      ProgName=Get_SuperName()
      If (ProgName(1:11).eq.'last_energy') DNG=.true.
*
*     Inside NUMERICAL_GRADIENT override input!
      If (ProgName(1:18).eq.'numerical_gradient') Then
         Call Put_iScalar('mp2prpt',0)
         DNG=.true.
         DoDens=.false.
         DoGrdt=.false.
      End If
*
      If(nSym .eq. 1) Then
         Call GetEnvF('EMIL_InLoop',emiloop)
         If (emiloop.eq.' ') emiloop='0'
         Call GetEnvF('MOLCAS_IN_GEO',inGeo)
         If ((emiloop(1:1).ne.'0') .and. inGeo(1:1) .ne. 'Y'
     &       .and. .not.DNG) Then
            Call Put_iScalar('mp2prpt',2)
            DoDens=.true.
            DoGrdt=.true.
         End If
      End If

      If (LumOrb) Then
         Lu_orb=7
         Call DaName(Lu_orb,'INPORB')
         l_Occup=nOrb(1)
         Do iSym=2,nSym
            l_Occup=l_Occup+nOrb(iSym)
         End Do
         Call GetMem('Occup','Allo','Real',ip_Occup,l_Occup)
         Call RDVEC('INPORB',Lu_orb,'COE',nSym,nBas,nOrb,CMO,
     &                       Work(ip_Occup),Eall,iDummy,
     &                       VecTitle,0,iErr)
         If (iErr.ne.0) Then
            Write(6,'(A,I4)') 'ERROR: RdVec returned code',iErr
            Call Abend()
         End If
         Call DaClos(Lu_orb)
         write(6,*)
         write(6,*) ' Input Orbitals read from INPORB: ',VecTitle
         If (.not.ERef_UsrDef) Then
            Write(6,*) ' WARNING: reference energy read from RunFile'
            Write(6,*) '          (may not correspond to orbitals)'
            Call Get_dScalar('SCF energy',Escf)
         End If
         write(6,*)
         iErr=0
         ip=ip_Occup-1
         Do iSym=1,nSym
            iCount=0
            Do i=1,nOrb(iSym)
               ip=ip+1
               If (abs(Work(ip)).gt.1.0d-14) iCount=iCount+1
            End Do
            If (iCount.ne.nOcc(iSym)) Then
               iErr=iErr+1
               Write(6,'(A,I2,A,I6,A,I6)')
     &         'WARNING: number of occupied orbitals in symmetry',iSym,
     &         ' is',iCount,' according to INPORB; from RunFile:',
     &         nOcc(iSym)
            End If
         End Do
         If (iErr.ne.0) Then
            Write(6,'(A,A)')
     &      'WARNING: occupation mismatch between RunFile and INPORB. ',
     &      'RunFile occupation will be used:'
            Write(6,'(8I6)') (nOcc(iSym),iSym=1,nSym)
            iErr=0
         End If
         Call GetMem('Occup','Free','Real',ip_Occup,l_Occup)
      Else
         Call Get_dScalar('SCF energy',Escf)
      EndIf
*
      If (lFre .and. nFre.gt.0) Then ! freeze lowest occupied orbitals
         FrePrt = iPrt.gt.0
         Call Freezer(Eall,nFre,nFro,nFro1,nOcc,nBas,nSym,FrePrt)
      End If
      Do iSym=1,nSym
         nFro(iSym)=nFro(iSym)+nFro1(iSym)
         nOcc(iSym)=nOcc(iSym)-nFro1(iSym)
         nDel(iSym)=nDel(iSym)+nDel1(iSym)
         nExt(iSym)=nExt(iSym)-nDel1(iSym)
         nOrb(iSym)=nOrb(iSym)-nDel1(iSym)
      End Do
      Call Put_iArray('nFroPT',NFRO,NSYM)
      Call Put_iArray('nDelPT',NDEL,NSYM)
      jOcc=0
      iExt=0
      iCount=0
      nOccT = 0
      nExtT = 0
      Do i= 1, nSym
         nOccT = nOccT + nOcc(i)
         nExtT = nExtT + nExt(i)
      End Do
*
      Do iSym=1,nSym
         iLow=iCount+nFro(iSym)+1
         iUpp=iLow+nOcc(iSym)-1
         Do i=iLow,iUpp
            jOcc=jOcc+1
            Eocc(jOcc)=Eall(i)
         End Do
         iLow=iUpp+1
         iUpp=iLow+nExt(iSym)-1
         Do i=iLow,iUpp
            iExt=iExt+1
            Eext(iExt)=Eall(i)
         End Do
         iCount=iCount+nBas(iSym)
      End Do
      iCount = 0
*
      If(DoDens) Then
         jFro = 0
         jDel = 0
         Do iSym = 1, nSym
*
            iLow = iCount+1
            iUpp = iLow+nFro(iSym)-1
            Do i = iLow, iUpp
               jFro=jFro+1
               EOcc(jFro+nOccT)= Eall(i)
            End Do
            iLow = iCount+nFro(iSym)+nOcc(iSym)+nExt(iSym)+1
            iUpp = iLow + nDel(iSym)-1
            Do i = iLow, iUpp
               jDel = jDel +1
               Eext(jDel+nExtT) = Eall(i)
            End Do
            iCount = iCount + nBas(iSym)
         End Do
      End If
      If ( lSFro .or. lSDel ) then
         LC=0
         LEO=0
         LEE=0
         LSQ=0
         Do iSym=1,nSym
            LC=LC+nBas(iSym)**2
            LSQ=Max(LSQ,nBas(iSym))
            LEO=LEO+nOcc(iSym)
            LEE=LEE+nExt(iSym)
         End Do
         Call GetMem('C','ALLO','REAL',IADC,LC)
         Call GetMem('EO','ALLO','REAL',IADEO,LEO)
         Call GetMem('EE','ALLO','REAL',IADEE,LEE)
         Call GetMem('SQ','ALLO','INTE',IADSQ,LSQ)
         call dcopy_(LC,CMO,1,WORK(IADC),1)
         call dcopy_(LEO,Eocc,1,WORK(IADEO),1)
         call dcopy_(LEE,Eext,1,WORK(IADEE),1)
         Call FrzDel(nFro2,iFro,Eocc,WORK(IADEO),
     &               nDel2,iDel,Eext,WORK(IADEE),
     &               CMO,WORK(IADC),IWORK(IADSQ))
         Call GetMem('SQ','FREE','INTE',IADSQ,LSQ)
         Call GetMem('EE','FREE','REAL',IADEE,LEE)
         Call GetMem('EO','FREE','REAL',IADEO,LEO)
         Call GetMem('LC','FREE','REAL',IADC,LC)
      End If
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
*
*     remove local copy of standard input
*
      Call Close_LuSpool(LuSpool)
*
      If (LumOrb .and. DoT1amp) Then
         Write (6,'(/,A)') 'ERROR!  Keywords incompatibility.'
         Write (6,'(/,A)') 'Both LUMOrb and T1AM were selected.'
         Write (6,'(/,A)') '***  I must shut down MBPT2 ! ***'
         Call Abend()
      EndIf
*
      If (all_Vir .and. DoMP2) Then
         Write (6,'(/,A)') 'WARNING!'
         Write (6,'(/,A)') 'Both VirAll and DoMP2 were selected.'
         Write (6,'(/,A)') '***  I turn off DoMP2 ! ***'
         DoMP2=.false.
      EndIf
*
*     turn off decomposition for parallel runs.
*
      If (DecoMP2 .and. ChoMP2_ChkPar()) Then
         Write (6,'(/,A)') 'WARNING!'
         Write (6,'(A)') 'Decomposition of MP2 integrals is not '
     &              //'possible. Turning off decomposition.'
         DecoMP2=.false.
      End If
*
      Call xFlush(6)
*
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
995   Write (6,*) 'RdInp: Error while reading input!'
      Write (6,'(A,A)') 'Last read line:',Line
      Call Abend()
      End
