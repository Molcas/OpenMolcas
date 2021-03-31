!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!***********************************************************************
      Subroutine Readinp_localisation()
!
!     Author: Y. Carissan [heavily modified by T.B. Pedersen].
!
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "debug.fh"
!
!TBP  Namelist /LOCALISATION/ dummy
!
      Character*20 SecNam
      Parameter (SecNam = 'Readinp_localisation')
!
      Character*180  Key, Line, Blank
      Character*180 Get_Ln
      External Get_Ln
!
      Integer  iPrintLevel
      External iPrintLevel
!
      Logical Thrs_UsrDef, LocModel_UsrDef, Freeze
      Parameter (ThrsDef = 1.0d-6, ThrRotDef = 1.0d-10)
      Parameter (ThrGradDef = 1.0d-2)
!
      LuSpool=17
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)
!
      Debug=.False.
!---- Locate "start of input"
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'LOCALISATION')
      Blank=' '
!
! Get print level
!
      iPL = iPrintLevel(-1)
!
! Default Parameters
!
      Do iSym = 1,nSym
         nOrb2Loc(iSym) = 0
         nFro(iSym) = 0
         nConstr(iSym) = 0
      End Do
      Skip = .False.
      LocOrb=Occupied
!     LocVir = .False.
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
!
! End Default Parameters
!
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
!
! DEBUg
!
 1000 Continue
      Debug=.True.
      Go To 999
!
! NORBitals
!
 2000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nOrb2Loc,nSym)
      nOrb2Loc_UsrDef=.True.
      Go To 999
!
! NFROzen
!
 2100 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I(1,nFro,nSym)
      nFro_UsrDef=.True.
      Freeze = .False.
      If (LocOrb.eq.Virtual) Then
!     If (LocVir) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'&
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '&
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
      EndIf
      Go To 999
!
! FREEze core orbitals (as defined by seward)
!
 2200 Continue
      If(.not.nFro_UsrDef)Then
        Call Get_iArray('Non valence orbitals',nFro,nSym)
        Freeze = .True.
        If (LocOrb.eq.Virtual) Then
!       If (LocVir) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'&
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '&
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
        EndIf
      End If
      Go To 999
!
! MAXImisation
!
 3000 Continue
      Maximisation=.True.
      Go To 999
!
! MINImisation
!
 3100 Continue
      Maximisation=.False.
      Go To 999
!
! NITErations or ITERations
!
 4000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_I1(1,NMxIter)
      Go To 999
!
! THREshold
!
 5000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,Thrs)
      Thrs_UsrDef = .True.
      Go To 999
!
! THRGrad
!
 5001 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ThrGrad)
      Go To 999
!
! CHOStart (use Cholesky orbitals as start guess for PM/Boys/ER)
!
 5100 Continue
      ChoStart = .True.
      Go To 999
!
! THRRotation
!
 6000 Continue
      Line=Get_Ln(LuSpool)
      Call Get_F1(1,ThrRot)
      Go To 999
!
! PIPEk-Mezey or PM
!
 7000 Continue
      LocModel = 1
      LocModel_UsrDef = .True.
      Go To 999
!
! BOYS
!
 7100 Continue
      LocModel = 2
      LocModel_UsrDef = .True.
      Go To 999
!
! CHOLesky
!
 7200 Continue
      LocModel = 3
      LocModel_UsrDef = .True.
      Go To 999
!
! EDMIston-Ruedenberg or ER
!
 7300 Continue
      LocModel = 4
      LocModel_UsrDef = .True.
      Go To 999
!
! SILEnt mode
!
 8000 Continue
      Silent = .True.
      Go To 999
!
! TEST localisation (orthonormality, density, etc.)
!
 9000 Continue
      Test_Localisation = .True.
      Go To 999
!
! ANALysis: generate bitmaps + histograms for density, original, and
!           local MOs. Analysis "per atom" is default, but does not work
!           with symmetry.
!
10000 Continue
      Analysis = .True.
      AnaAtom  = nSym.eq.1
      Go To 999
!
! ANAAtom: generate bitmaps + histograms for density, original, and
!          local MOs. Analysis "per atom".
!
10100 Continue
      Analysis = .True.
      AnaAtom  = .True.
      Go To 999
!
! ANAShell: generate bitmaps + histograms for density, original, and
!           local MOs. Analysis "per shell".
!
10200 Continue
      Analysis = .True.
      AnaAtom  = .False.
      Go To 999
!
! MAX norm: use max element as norm in analysis.
!
10300 Continue
      AnaNrm  = 'Max'
      Go To 999
!
! FROBenius norm: use Frobenius norm in analysis.
!
10310 Continue
      AnaNrm  = 'Fro'
      Go To 999
!
! NOMO: do not print MOs.
!
11000 Continue
      PrintMOs = .False.
      Go To 999
!
! TIME: time the localisation procedure explicitly
!
12000 Continue
      Timing = .True.
      Go To 999
!
! NOTI: do not time the localisation procedure explicitly
!
12100 Continue
      Timing = .False.
      Go To 999
!
! ERFUnctional: evaluate Edmiston-Ruedenberg functional with original
!               and localised orbitals.
!
13000 Continue
      EvalER = .True.
      Go To 999
!
! ORDEr: order localised orbitals according to Cholesky orbital
!        ordering, defined according to the overlap between the two sets
!        orbitals.
!
14000 Continue
      Order = .True.
      Go To 999
!
! VIRTual: localise virtual orbitals
!
15000 Continue
      LocOrb = Virtual
!     LocVir = .True.
      If (nFro_UsrDef .or. Freeze) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'You have chosen to freeze some orbitals AND asked'&
     &             //' to localize the virtual space. This implies that'
         Write (6,*) ' some virtual orbitals will be kept frozen, thus '&
     &             //'they will be left unchanged. Is that OK ?'
         Write (6,*)
      EndIf
      Go To 999
!
! OCCUpied: localise occupied orbitals
!
15100 Continue
      LocOrb = Occupied
!     LocVir = .False.
      Go To 999
!
! ALL: localise all orbitals
!
15200 Continue
      LocOrb = All
      Go To 999
!
! PAO : compute projected AOs that span the virtual space
!       (= all - occupied - frozen) using Cholesky decomposition to
!       remove linear dependence.
!
16000 Continue
      LocPAO = .True.
      LocModel = 3
      LocModel_UsrDef = .True.
      Go To 999
!
! ANAP: Special analysis for Cholesky PAOs before orthonormalization.
!
16100 Continue
      AnaPAO = .True.
      Go To 999
!
! DOMAin: set up orbital domains (Pulay-style) and pair domains
!
17000 Continue
      DoDomain = .True.
      Go To 999
!
! THRDomain: thresholds for setting up orbital domains.
!            First value is the minimum sum of gross atomic Mulliken
!            charge for each orbital.
!            Second value is the threshold used for the completeness
!            check of Boughton and Pulay.
!
17100 Continue
      DoDomain = .True.
      Line=Get_Ln(LuSpool)
      Call Get_F(1,ThrDomain,2)
      Go To 999
!
! THRPairdomain: thresholds for setting up pair domains.
!                3 values needed: Rs, Rw, Rd (in bohr).
!                Let R be the minimum distance between any two atoms in
!                a pair domain. Then,
!                     R <= Rs: strong pair
!                Rs < R <= Rw: weak pair
!                Rw < R <= Rd: distant pair
!                Rd < R <= Rd: very distant pair
!
17200 Continue
      DoDomain = .True.
      Line=Get_Ln(LuSpool)
      Call Get_F(1,ThrPairDomain,3)
      Go To 999
!
! ANADomain: analysis of orbital domains
!
17300 Continue
      AnaDomain = .True.
      DoDomain = .True.
      Go To 999
!
! SKIP: skip localisation completely; i.e. only perform analysis etc.
!
18000 Continue
      Skip = .True.
      Go To 999
!
! LOCN: localized natural orbitals
!
18100 Continue
      LocNatOrb=.true.
      Read(LuSpool,*) nActa, ThrSel
!     Read atom names
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
!
! LOCC: localized canonical orbitals
!
18200 Continue
      LocCanOrb=.true.
      Read(LuSpool,*) nActa, ThrSel
!     Read atom names
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
!
! WAVE: wavelet transform of the MO basis
!
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
!
! CONS: constrained natural orbital analysis
!
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
!
! FILE: filename with input orbitals
!
18500 Continue
      Line=Get_Ln(LuSpool)
!     Filename is read in localisation.f
      Go To 999
!
! END of Input
!
99999 Continue

!     ==============
!     Postprocessing
!     ==============

!     Only Cholesky localisation can be run with symmetry.
!     Reset if the user explicitly requested a localisation
!     model different from Cholesky.
!     ------------------------------------------------------------------

      If (nSym.gt.1 .and. LocModel_UsrDef .and. .not.LocModel.eq.3) Then
         Write (6,*)
         Write (6,*) 'WARNING!!!'
         Write (6,*) 'The localisation model you have suggested in '    &
     &             //'combination with symmetry will not work.'
         Write (6,*) 'In this case the program defaults to the '        &
     &             //'Cholesky localisation scheme!'
         Write (6,*)
         LocModel = 3
      End If

!     Special settings for PAO: LocModel must be 3 (i.e. Cholesky), else
!     we cancel PAO. For PAO, special analysis may be activated by the
!     DEBUg keyword.
!     ------------------------------------------------------------------

      LocPAO = LocPAO .and. LocModel.eq.3
      If (LocPAO) Then
         AnaPAO = AnaPAO .or. Debug
      Else
         AnaPAO = .False.
      End If

!     Set orbitals to localize.
!     If the user specified which orbitals to localise, frozen orbitals
!     is not our problem (the user must set those as well or use the
!     default).
!     If the user did not specify which orbitals to localise, we have
!     several cases (remembering that virtual localisation is default
!     for the PAO method):
!     For virtual localisation (specified with keyword VIRTual):
!        user-defined frozen orbitals are regarded as frozen virtual,
!        else we freeze all occupied and localise the rest. Note that
!        keyword FREEze is ignored for virtual localisation!
!     For occupied localisation (default for all methods but PAO):
!        user-defined frozen orbitals (possibly by FREEze keyword)
!        are subtracted from the set of occupied orbitals, else we
!        localise all occupied.
!     -----------------------------------------------------------------

      If (nOrb2Loc_UsrDef) Then
         LocOrb = Occupied
!        LocVir = .False.
      Else
!        If (LocOrb.eq.Virtual .or. LocPAO) Then
!           LocOrb = Virtual
!        End If
!        LocVir = LocVir .or. LocPAO ! virtual is default for PAO
         If (LocOrb.eq.Virtual) Then ! virtual localisation
!        If (LocVir) Then ! virtual localisation
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
!        Else ! occupied localisation
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
                  nOrb2Loc(iSym) = nOccInp(iSym) + nVirInp(iSym)        &
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

!     If Cholesky, reset default threshold (unless user defined).
!     -----------------------------------------------------------

      If (LocModel.eq.3 .and. .not.Thrs_UsrDef) Then
         Thrs = 1.0d-8
      End If

!     No need to order Cholesky MOs.
!     ------------------------------

      Order = Order .and. LocModel.ne.3

!     No need to evaluate ER function for Edmiston-Ruedenberg.
!     --------------------------------------------------------

      EvalER = EvalER .and. LocModel.ne.4

!     Turn on localisation test if debug or analysis was specified.
!     -------------------------------------------------------------

      Test_Localisation = Test_Localisation .or. Debug .or. Analysis

      Go To 9950

9940  Continue
      WRITE(6,*)' READIN: Premature end of file when reading selected'
      WRITE(6,*)' atoms in keyword LOCN'
      CALL ABEND()

9950  Continue

      End
