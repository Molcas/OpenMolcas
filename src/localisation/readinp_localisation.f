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
* Copyright (C) Yannick Carissan                                       *
*               Thomas Bondo Pedersen                                  *
************************************************************************
      Subroutine Readinp_localisation()
c
c     Author: Y. Carissan [heavily modified by T.B. Pedersen].
c
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "debug.fh"
c
cTBP  Namelist /LOCALISATION/ dummy
c
      Character*20 SecNam
      Parameter (SecNam = 'Readinp_localisation')
c
      Character*180  Key, Line, Blank
      Character*180 Get_Ln
      External Get_Ln
c
      Integer  iPrintLevel
      External iPrintLevel
c
      Logical Thrs_UsrDef, LocModel_UsrDef, Freeze
      Parameter (ThrsDef = 1.0d-6, ThrRotDef = 1.0d-10)
      Parameter (ThrGradDef = 1.0d-2)
c
      LuSpool=17
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)
c
      Debug=.False.
c---- Locate "start of input"
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'LOCALISATION')
      Blank=' '
c
c Get print level
c
      iPL = iPrintLevel(-1)
c
c Default Parameters
c
      Do iSym = 1,nSym
         nOrb2Loc(iSym) = 0
         nFro(iSym) = 0
         nConstr(iSym) = 0
      End Do
      Skip = .False.
      LocOrb=Occupied
*     LocVir = .False.
      Thrs_UsrDef = .False.
      nOrb2Loc_UsrDef = .False.
      nFro_UsrDef = .False.
      Freeze = .False.
      If (iPL .ge. 4) Then
         Debug = .True.
      Else
         Debug = .False.
      End If
      Maximisation = .True.
      ChoStart = .False.
      If (iPL .lt. 3) Then
         Silent = .True.
      Else
         Silent = .False.
      End If
      LocModel = 1  ! Pipek-Mezey localisation
      If (nSym.gt.1) LocModel = 3  ! Cholesky localisation
      LocModel_UsrDef = .False.
      Test_Localisation = .False.
      NMxIter = 300
      Thrs = ThrsDef
      ThrRot = ThrRotDef
      ThrGrad = ThrGradDef
      Analysis = .False.
      AnaAtom  = nSym.eq.1
      AnaNrm = 'Fro'
      PrintMOs = .True.
      Timing = .True.
      EvalER = .False.
      Order = .False.
      LocPAO = .False.
      AnaPAO = .False.
      AnaPAO_Save = AnaPAO
      DoDomain = .False.
      AnaDomain = .False.
      ThrDomain(1) = 9.0d-1
      ThrDomain(2) = 2.0d-2
      ThrPairDomain(1) = 1.0d-10
      ThrPairDomain(2) = 1.0d1
      ThrPairDomain(3) = 1.5d1
      LocNatOrb=.false.
      LocCanOrb=.false.
      Wave=.false.
      iWave=0
      DoCNOs=.false.
c
c End Default Parameters
c
  999 Continue
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)
      If (Line(1:4).eq.'DEBU') Go To 1000
      If (Line(1:4).eq.'NORB') Go To 2000
      If (Line(1:4).eq.'NFRO') Go To 2100
      If (Line(1:4).eq.'FREE') Go To 2200
      If (Line(1:4).eq.'MAXI') Go To 3000
      If (Line(1:4).eq.'MINI') Go To 3100
      If (Line(1:4).eq.'NITE') Go To 4000
      If (Line(1:4).eq.'ITER') Go To 4000
      If (Line(1:4).eq.'THRE') Go To 5000
      If (Line(1:4).eq.'THRG') Go To 5001
      If (Line(1:4).eq.'CHOS') Go To 5100
      If (Line(1:4).eq.'THRR') Go To 6000
      If (Line(1:4).eq.'PIPE') Go To 7000
      If (Line(1:4).eq.'PM  ') Go To 7000
      If (Line(1:4).eq.'BOYS') Go To 7100
      If (Line(1:4).eq.'CHOL') Go To 7200
      If (Line(1:4).eq.'EDMI') Go To 7300
      If (Line(1:4).eq.'ER  ') Go To 7300
      If (Line(1:4).eq.'SILE') Go To 8000
      If (Line(1:4).eq.'TEST') Go To 9000
      If (Line(1:4).eq.'ANAL') Go To 10000
      If (Line(1:4).eq.'ANAA') Go To 10100
      If (Line(1:4).eq.'ANAS') Go To 10200
      If (Line(1:4).eq.'MAX ') Go To 10300
      If (Line(1:4).eq.'FROB') Go To 10310
      If (Line(1:4).eq.'NOMO') Go To 11000
      If (Line(1:4).eq.'TIME') Go To 12000
      If (Line(1:4).eq.'NOTI') Go To 12100
      If (Line(1:4).eq.'ERFU') Go To 13000
      If (Line(1:4).eq.'ORDE') Go To 14000
      If (Line(1:4).eq.'VIRT') Go To 15000
      If (Line(1:4).eq.'OCCU') Go To 15100
      If (Line(1:4).eq.'ALL ') Go To 15200
      If (Line(1:4).eq.'PAO ') Go To 16000
      If (Line(1:4).eq.'ANAP') Go To 16100
      If (Line(1:4).eq.'DOMA') Go To 17000
      If (Line(1:4).eq.'THRD') Go To 17100
      If (Line(1:4).eq.'THRP') Go To 17200
      If (Line(1:4).eq.'ANAD') Go To 17300
      If (Line(1:4).eq.'SKIP') Go To 18000
      If (Line(1:4).eq.'LOCN') Go To 18100
      If (Line(1:4).eq.'LOCC') Go To 18200
      If (Line(1:4).eq.'WAVE') Go To 18300
      If (Line(1:4).eq.'CONS') Go To 18400
      If (Line(1:4).eq.'FILE') Go To 18500
      If (Line(1:4).eq.'END ') Go To 99999
      Write (6,*) 'Unidentified key word  : ',Key
      Write (6,*) 'Internal representation: ',Line(1:4)
      Call FindErrorLine
      Call Quit_OnUserError()
c
c DEBUg
c
 1000 Continue
      Debug=.True.
      Go To 999
c
c NORBitals
c
 2000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nOrb2Loc,nSym)
      nOrb2Loc_UsrDef=.True.
      Go To 999
c
c NFROzen
c
 2100 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nFro,nSym)
      nFro_UsrDef=.True.
      Freeze = .False.
      If (LocOrb.eq.Virtual) Then
*     If (LocVir) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
      EndIf
      Go To 999
c
c FREEze core orbitals (as defined by seward)
c
 2200 Continue
      If(.not.nFro_UsrDef)Then
        Call Get_iArray('Non valence orbitals',nFro,nSym)
        Freeze = .True.
        If (LocOrb.eq.Virtual) Then
*       If (LocVir) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
        EndIf
      End If
      Go To 999
c
c MAXImisation
c
 3000 Continue
      Maximisation=.True.
      Go To 999
c
c MINImisation
c
 3100 Continue
      Maximisation=.False.
      Go To 999
c
c NITErations or ITERations
c
 4000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,NMxIter)
      Go To 999
c
c THREshold
c
 5000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,Thrs)
      Thrs_UsrDef = .True.
      Go To 999
c
c THRGrad
c
 5001 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ThrGrad)
      Go To 999
c
c CHOStart (use Cholesky orbitals as start guess for PM/Boys/ER)
c
 5100 Continue
      ChoStart = .True.
      Go To 999
c
c THRRotation
c
 6000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ThrRot)
      Go To 999
c
c PIPEk-Mezey or PM
c
 7000 Continue
      LocModel = 1
      LocModel_UsrDef = .True.
      Go To 999
c
c BOYS
c
 7100 Continue
      LocModel = 2
      LocModel_UsrDef = .True.
      Go To 999
c
c CHOLesky
c
 7200 Continue
      LocModel = 3
      LocModel_UsrDef = .True.
      Go To 999
c
c EDMIston-Ruedenberg or ER
c
 7300 Continue
      LocModel = 4
      LocModel_UsrDef = .True.
      Go To 999
c
c SILEnt mode
c
 8000 Continue
      Silent = .True.
      Go To 999
c
c TEST localisation (orthonormality, density, etc.)
c
 9000 Continue
      Test_Localisation = .True.
      Go To 999
c
c ANALysis: generate bitmaps + histograms for density, original, and
c           local MOs. Analysis "per atom" is default, but does not work
c           with symmetry.
c
10000 Continue
      Analysis = .True.
      AnaAtom  = nSym.eq.1
      Go To 999
c
c ANAAtom: generate bitmaps + histograms for density, original, and
c          local MOs. Analysis "per atom".
c
10100 Continue
      Analysis = .True.
      AnaAtom  = .True.
      Go To 999
c
c ANAShell: generate bitmaps + histograms for density, original, and
c           local MOs. Analysis "per shell".
c
10200 Continue
      Analysis = .True.
      AnaAtom  = .False.
      Go To 999
c
c MAX norm: use max element as norm in analysis.
c
10300 Continue
      AnaNrm  = 'Max'
      Go To 999
c
c FROBenius norm: use Frobenius norm in analysis.
c
10310 Continue
      AnaNrm  = 'Fro'
      Go To 999
c
c NOMO: do not print MOs.
c
11000 Continue
      PrintMOs = .False.
      Go To 999
c
c TIME: time the localisation procedure explicitly
c
12000 Continue
      Timing = .True.
      Go To 999
c
c NOTI: do not time the localisation procedure explicitly
c
12100 Continue
      Timing = .False.
      Go To 999
c
c ERFUnctional: evaluate Edmiston-Ruedenberg functional with original
c               and localised orbitals.
c
13000 Continue
      EvalER = .True.
      Go To 999
c
c ORDEr: order localised orbitals according to Cholesky orbital
c        ordering, defined according to the overlap between the two sets
c        orbitals.
c
14000 Continue
      Order = .True.
      Go To 999
c
c VIRTual: localise virtual orbitals
c
15000 Continue
      LocOrb = Virtual
*     LocVir = .True.
      If (nFro_UsrDef .or. Freeze) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
      EndIf
      Go To 999
c
c OCCUpied: localise occupied orbitals
c
15100 Continue
      LocOrb = Occupied
*     LocVir = .False.
      Go To 999
c
c ALL: localise all orbitals
c
15200 Continue
      LocOrb = All
      Go To 999
c
c PAO : compute projected AOs that span the virtual space
c       (= all - occupied - frozen) using Cholesky decomposition to
c       remove linear dependence.
c
16000 Continue
      LocPAO = .True.
      LocModel = 3
      LocModel_UsrDef = .True.
      Go To 999
c
c ANAP: Special analysis for Cholesky PAOs before orthonormalization.
c
16100 Continue
      AnaPAO = .True.
      Go To 999
c
c DOMAin: set up orbital domains (Pulay-style) and pair domains
c
17000 Continue
      DoDomain = .True.
      Go To 999
c
c THRDomain: thresholds for setting up orbital domains.
c            First value is the minimum sum of gross atomic Mulliken
c            charge for each orbital.
c            Second value is the threshold used for the completeness
c            check of Boughton and Pulay.
c
17100 Continue
      DoDomain = .True.
      Line=Get_Ln(LuSpool)
      Call Get_F(1,ThrDomain,2)
      Go To 999
c
c THRPairdomain: thresholds for setting up pair domains.
c                3 values needed: Rs, Rw, Rd (in bohr).
c                Let R be the minimum distance between any two atoms in
c                a pair domain. Then,
c                     R <= Rs: strong pair
c                Rs < R <= Rw: weak pair
c                Rw < R <= Rd: distant pair
c                Rd < R <= Rd: very distant pair
c
17200 Continue
      DoDomain = .True.
      Line=Get_Ln(LuSpool)
      Call Get_F(1,ThrPairDomain,3)
      Go To 999
c
c ANADomain: analysis of orbital domains
c
17300 Continue
      AnaDomain = .True.
      DoDomain = .True.
      Go To 999
c
c SKIP: skip localisation completely; i.e. only perform analysis etc.
c
18000 Continue
      Skip = .True.
      Go To 999
c
c LOCN: localized natural orbitals
c
18100 Continue
      LocNatOrb=.true.
      Read(LuSpool,*) nActa, ThrSel
*     Read atom names
18101 Read(LuSpool,'(A)',End=9940) Line
      If ( Line(1:1).eq.'*' ) Goto 18101
      If ( Line.eq.Blank ) Goto 18101
      Call UpCase(Line)
      Do i=1,nActa
       Call LeftAd(Line)
       If( Line.eq.Blank ) Goto 9940
       j=Index(Line,' ')
       NamAct(i)=Line(1:j-1)
       Line(1:j-1)=Blank(1:j-1)
      Enddo
      Go To 999
c
c LOCC: localized canonical orbitals
c
18200 Continue
      LocCanOrb=.true.
      Read(LuSpool,*) nActa, ThrSel
*     Read atom names
      Read(LuSpool,'(A)',End=9940) Line
      If ( Line(1:1).eq.'*' ) Goto 18101
      If ( Line.eq.Blank ) Goto 18101
      Call UpCase(Line)
      Do i=1,nActa
       Call LeftAd(Line)
       If( Line.eq.Blank ) Goto 9940
       j=Index(Line,' ')
       NamAct(i)=Line(1:j-1)
       Line(1:j-1)=Blank(1:j-1)
      Enddo
      Go To 999
c
c WAVE: wavelet transform of the MO basis
c
18300 Continue
      Read(LuSpool,*) iWave
      If (iWave.ne.0 .and. iWave.ne.1) Then
         Write (6,*) ' WARNING!!!'
         Write (6,*) ' Incorrect bit-switch input parameter for WAVE.'
         Write (6,*) ' I will continue with the default value.'
         iWave=0
      EndIf
      Skip=.false.
      Wave=.true.
      LocModel=0
      Go To 999
c
c CONS: constrained natural orbital analysis
c
18400 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nConstr,nSym)
      MxConstr=0
      Do i=1,nSym
         MxConstr=Max(MxConstr,nConstr(i))
      End Do
      If (MxConstr.gt.16) Then
        Write (6,*) ' ERROR  !!!'
        Write (6,*) ' Cannot handle more than 16 constraints/symm !!'
        Write (6,*) ' Increase size of indxC in Get_CNOs and recompile.'
      EndIf
      DoCNOs=.true.
      Skip=.false.
      PrintMOs=.false.
      LocModel=0
      Go To 999
c
c FILE: filename with input orbitals
c
18500 Continue
      Line=Get_Ln(LuSpool)
C     Filename is read in localisation.f
      Go To 999
c
c END of Input
c
99999 Continue

C     ==============
C     Postprocessing
C     ==============

C     Only Cholesky localisation can be run with symmetry.
C     Reset if the user explicitly requested a localisation
C     model different from Cholesky.
C     ------------------------------------------------------------------

      If (nSym.gt.1 .and. LocModel_UsrDef .and. .not.LocModel.eq.3) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'The localisation model you have suggested in '
     &             //'combination with symmetry will not work.'
         Write (6,*) 'In this case the program defaults to the '
     &             //'Cholesky localisation scheme!'
         Write (6,*)
         LocModel = 3
      End If

C     Special settings for PAO: LocModel must be 3 (i.e. Cholesky), else
C     we cancel PAO. For PAO, special analysis may be activated by the
C     DEBUg keyword.
C     ------------------------------------------------------------------

      LocPAO = LocPAO .and. LocModel.eq.3
      If (LocPAO) Then
         AnaPAO = AnaPAO .or. Debug
      Else
         AnaPAO = .False.
      End If

C     Set orbitals to localize.
C     If the user specified which orbitals to localise, frozen orbitals
C     is not our problem (the user must set those as well or use the
C     default).
C     If the user did not specify which orbitals to localise, we have
C     several cases (remembering that virtual localisation is default
C     for the PAO method):
C     For virtual localisation (specified with keyword VIRTual):
C        user-defined frozen orbitals are regarded as frozen virtual,
C        else we freeze all occupied and localise the rest. Note that
C        keyword FREEze is ignored for virtual localisation!
C     For occupied localisation (default for all methods but PAO):
C        user-defined frozen orbitals (possibly by FREEze keyword)
C        are subtracted from the set of occupied orbitals, else we
C        localise all occupied.
C     -----------------------------------------------------------------

      If (nOrb2Loc_UsrDef) Then
         LocOrb = Occupied
*        LocVir = .False.
      Else
*        If (LocOrb.eq.Virtual .or. LocPAO) Then
*           LocOrb = Virtual
*        End If
*        LocVir = LocVir .or. LocPAO ! virtual is default for PAO
         If (LocOrb.eq.Virtual) Then ! virtual localisation
*        If (LocVir) Then ! virtual localisation
            If (nFro_UsrDef) Then
               Do iSym = 1,nSym
                  nFro(iSym) = nOccInp(iSym) + nFro(iSym)
                  nOrb2Loc(iSym) = nOrb(iSym) - nFro(iSym)
               End Do
            Else
               Do iSym = 1,nSym
                  nFro(iSym) = nOccInp(iSym)
                  nOrb2Loc(iSym) = nVirInp(iSym)
               End Do
            End If
         Else If (LocOrb.eq.Occupied) Then ! occupied localisation
*        Else ! occupied localisation
            If (nFro_UsrDef .or. Freeze) Then
               Do iSym = 1,nSym
                  nOrb2Loc(iSym) = nOccInp(iSym) - nFro(iSym)
               End Do
            Else
               Do iSym = 1,nSym
                  nOrb2Loc(iSym) = nOccInp(iSym)
               End Do
            End If
         Else ! All
            If (nFro_UsrDef .or. Freeze) Then
               Do iSym = 1,nSym
                  nOrb2Loc(iSym) = nOccInp(iSym) + nVirInp(iSym)
     &                           - nFro(iSym)
               End Do
            Else
               Do iSym = 1,nSym
                  nOrb2Loc(iSym) = nOccInp(iSym) + nVirInp(iSym)
               End Do
            End If
         End If
      End If

      If (Debug) Then
         Write(6,'(/,A,A)') SecNam,': orbital definitions:'
         Write(6,'(A,8I9)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
         Write(6,'(A,8I9)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
         Write(6,'(A,8I9)') 'nOccInp : ',(nOccInp(iSym),iSym=1,nSym)
         Write(6,'(A,8I9)') 'nVirInp : ',(nVirInp(iSym),iSym=1,nSym)
         Write(6,'(A,8I9)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
         Write(6,'(A,8I9,/)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
      End If

C     If Cholesky, reset default threshold (unless user defined).
C     -----------------------------------------------------------

      If (LocModel.eq.3 .and. .not.Thrs_UsrDef) Then
         Thrs = 1.0d-8
      End If

C     No need to order Cholesky MOs.
C     ------------------------------

      Order = Order .and. LocModel.ne.3

C     No need to evaluate ER function for Edmiston-Ruedenberg.
C     --------------------------------------------------------

      EvalER = EvalER .and. LocModel.ne.4

C     Turn on localisation test if debug or analysis was specified.
C     -------------------------------------------------------------

      Test_Localisation = Test_Localisation .or. Debug .or. Analysis

      Go To 9950

9940  Continue
      WRITE(6,*)' READIN: Premature end of file when reading selected'
      WRITE(6,*)' atoms in keyword LOCN'
      CALL ABEND()

9950  Continue

      End
