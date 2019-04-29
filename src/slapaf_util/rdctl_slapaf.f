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
      Subroutine RdCtl_Slapaf(iRow,iInt,nFix,LuSpool,Dummy_Call)
      use AI
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "nadc.fh"
#include "weighting.fh"
#include "print.fh"
      Logical Found, Dummy_Call
      Character*8 Command
      Character*180 Get_Ln
      Character*16 FilNam
      Character*3 MEPLab
      External Get_SuperName
      Character*100 Get_SuperName
      Character*100 SuperName
*
*     Compare with inputil.f. Note that here Line is defined in
*     info:slapaf.fh. Otherwise the common cgetln should be
*     identical in size.
*
* mxn should be len(line)/2+1
      parameter (mxn=91)
      common/cgetln/ ncol, jstrt(mxn),jend(mxn)
      Integer StrnLn
      External Get_Ln, StrnLn
      Logical External_UDC, External_Case,
     &        Explicit_IRC, Expert, ThrInp, FirstNum
#include "angstr.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      Call QEnter('RdCtl_Slapaf')
      Expert=.False.
      Lu=6
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Initiate some parameters
*
      Call Init_Slapaf(iRow)
      iPrint=nPrint(iRout)
      iSetAll=2**30 - 1
*
      Call f_Inquire('UDC.Gateway',External_UDC)
      LuRd2=LuSpool
      External_Case=.False.
*
      iMEP=0
      Explicit_IRC=.False.
      lCtoF=.False.
      WeightedConstraints=.False.
      ThrInp=.False.
      Call Qpg_iScalar('nMEP',Found)
      If (Found) Call Get_iScalar('nMEP',iMEP)
      If (iMEP.eq.0) Then
         iOff_Iter=0
         Call Put_iScalar('iOff_Iter',iOff_Iter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     When called from outside Slapaf or as a dummy call, process no
*     input but proceed with default values only.
*
      SuperName= Get_Supername()
      If ((SuperName.ne.'slapaf').or.Dummy_Call) Then
         Char='END '
         Go To 666
      End If
*                                                                      *
************************************************************************
**************************   Input section   ***************************
************************************************************************
*                                                                      *
      LuRd=LuSpool
      Call RdNlst(LuRd,'SLAPAF')
      Command='&SLAPAF'
 999  Char=Get_Ln(LuRd)
 666  Continue
      Call UpCase(Char)
      Command=Char(1:8)
C     Write (Lu,'(A)') Char
C     Write (Lu,*) iOptC
      If (Char.eq.BLine) Go To 999
      If (Char(1:1).eq.'*') Go To 999
      If (Char(1:4).eq.'AI  ') Go To 100
      If (Char(1:4).eq.'AIAM') Go To 101
      If (Char(1:4).eq.'AIL ') Go To 102
      If (Char(1:4).eq.'AINX') Go To 103
      If (Char(1:4).eq.'AIP ') Go To 104
      If (Char(1:4).eq.'AISP') Go To 105
      If (Char(1:4).eq.'AIMI') Go To 106
      If (Char(1:4).eq.'AIME') Go To 107
      If (Char(1:4).eq.'BAKE') Go To 926
      If (Char(1:4).eq.'C1-D') Go To 936
      If (Char(1:4).eq.'C2-D') Go To 937
      If (Char(1:4).eq.'CART') Go To 918
      If (Char(1:4).eq.'CNWE') Go To 990
      If (Char(1:4).eq.'CONS') Go To 9478
      If (Char(1:4).eq.'CTOF') Go To 904
      If (Char(1:4).eq.'CUBI') Go To 947
      If (Char(1:4).eq.'DDVS') Go To 9271
      If (Char(1:4).eq.'DELT') Go To 946
      If (Char(1:4).eq.'DISO') Go To 9452
      If (Char(1:4).eq.'DXDX') Go To 939
      If (Char(1:4).eq.'DXG ') Go To 940
      If (Char(1:4).eq.'END ') Go To 998
      If (Char(1:4).eq.'EXPE') Go To 993
      If (Char(1:4).eq.'EXTR') Go To 971
      If (Char(1:4).eq.'FALC') Go To 800
      If (Char(1:4).eq.'FIND') Go To 963
      If (Char(1:4).eq.'FUZZ') Go To 123
      If (Char(1:4).eq.'GDX ') Go To 940
      If (Char(1:4).eq.'GG  ') Go To 941
      If (Char(1:4).eq.'GNRM') Go To 968
      If (Char(1:4).eq.'GRAD') Go To 979
      If (Char(1:4).eq.'HRMS') Go To 995
      If (Char(1:4).eq.'HUPD') Go To 914
      If (Char(1:4).eq.'HWRS') Go To 929
      If (Char(1:4).eq.'INTE') Go To 902
      If (Char(1:4).eq.'IRC ') Go To 997
      If (Char(1:4).eq.'ITER') Go To 925
      If (Char(1:4).eq.'LAST') Go To 9280
      If (Char(1:4).eq.'LINE') Go To 9281
      If (Char(1:4).eq.'MAXS') Go To 915
      If (Char(1:4).eq.'MEP-'.or. Char(1:4).eq.'MEP ') Go To 964
      If (Char(1:4).eq.'MEPA'.or. Char(1:4).eq.'IRCA') Go To 322
      If (Char(1:4).eq.'MEPS'.or. Char(1:4).eq.'IRCS') Go To 9971
      If (Char(1:4).eq.'MEPT'.or. Char(1:4).eq.'IRCT') Go To 321
      If (Char(1:4).eq.'MODE') Go To 942
      If (Char(1:4).eq.'NMEP'.or. Char(1:4).eq.'NIRC') Go To 965
      If (Char(1:4).eq.'NEWT') Go To 935
      If (Char(1:4).eq.'NOEM') Go To 991
      If (Char(1:4).eq.'NOHW') Go To 960
      If (Char(1:4).eq.'NOLA') Go To 930
      If (Char(1:4).eq.'NOLI') Go To 928
      If (Char(1:4).eq.'NOWB') Go To 984
      If (Char(1:4).eq.'NOWC') Go To 985
      If (Char(1:4).eq.'NOWH') Go To 986
      If (Char(1:4).eq.'NUME') Go To 945
      If (Char(1:4).eq.'OLDF') Go To 903
      If (Char(1:4).eq.'PRFC') Go To 9201
      If (Char(1:4).eq.'PRIN') Go To 920
      If (Char(1:4).eq.'RATI') Go To 938
      If (Char(1:4).eq.'REDU') Go To 994
      If (Char(1:4).eq.'REAC') Go To 996
      If (Char(1:4).eq.'REFE') Go To 966
      If (Char(1:4).eq.'RHID') Go To 988
      If (Char(1:4).eq.'RMEP') Go To 980
      If (Char(1:4).eq.'RS-P') Go To 967
      If (Char(1:4).eq.'RTRN') Go To 962
      If (Char(1:4).eq.'SCHL') Go To 927
      If (Char(1:4).eq.'SUPS') Go To 911
      If (Char(1:4).eq.'THER') Go To 9451
      If (Char(1:4).eq.'THRS') Go To 908
      If (Char(1:4).eq.'TOLE') Go To 909
      If (Char(1:4).eq.'TRAC') Go To 910
      If (Char(1:4).eq.'TS  ') Go To 951
      If (Char(1:4).eq.'TSCO') Go To 320
      If (Char(1:4).eq.'VDWB') Go To 981
      If (Char(1:4).eq.'VDWC') Go To 982
      If (Char(1:4).eq.'VDWH') Go To 983
      If (Char(1:4).eq.'WIND') Go To 934
      Call WarningMessage(2,'Error in RdCtl_Slapaf')
      If (Char(1:1).eq.' ') Then
         Write (Lu,*) ' RdCtl_Slapaf: Command line starts with a blank.'
      Else
         Write (Lu,*)
         Write (Lu,*) ' *********** ERROR ***********'
         Write (Lu,*) ' The program has been supplied'
         Write (Lu,*) ' with an unknown command.     '
         Write (Lu,*) ' *****************************'
      End If
      Write (Lu,'(A)') Char
      Call Quit_OnUserError()
*                                                                      *
****** INTE ************************************************************
*                                                                      *
*     Read the internal coordinate specification.
*
 902  Continue
      New_Line=1
      Lu_UDIC=91
      FilNam='UDIC'
cc      Open(Lu_UDIC,File=FilNam,Form='FORMATTED',Status='UNKNOWN')
      call molcas_open(Lu_UDIC,FilNam)
      ReWind(Lu_UDIC)
*
*     mInt is the number of internal coordinates you will define.
*     Subroutine DefInt defines the B matrix.
*     The matrix B relates a shift in an internal coordinate to
*     shifts in cartesian coordinates,
*
*               |dq> = B |dx>
*                      =
*     and has the dimension (3*nsAtom x mInt).
 992  Continue
         Line=Get_Ln(LuRd)
         Call UpCase(Line)
         If (Line(1:4).eq.'END ') Then
            Close(Lu_UDIC)
            Go To 999
         End If
*
*        Here is a fix because auto will break up the lines if there is an
*        equal sign in the input.
*
*        Lines with VARY or FIX doesn't have equal signs
*
         If (Line(1:4).eq.'VARY') nBVec=iRow
         If (Line(1:4).eq.'VARY'.or.
     &       Line(1:3).eq.'FIX' .or.
     &       Line(1:4).eq.'ROWH') Then
            New_Line=0
         End If
*
 111     Continue
         If (New_Line.eq.1) Then
            If (Index(Line,'=').eq.0) Call FixEqualSign2(Line,LuRd,
     &                                                   Lu_UDIC,iRow,
     &                                                   New_Line)
            If (New_Line.eq.2) Then
               Close(Lu_UDIC)
               Go To 999
            End If
            Go To 111
         End If
*
         iRow = iRow + 1
*
         Write (Lu_UDIC,'(A)') Line
*
*        If this line does not have a continuation the next line should
*        have a equal sign!
         If (Index(Line,'&').eq.0) New_Line=1
      Go To 992
*                                                                      *
****** CTOF ************************************************************
*                                                                      *
*     Read the internal (C)oordinate specification (TO) be (F)ollowed.
*
 904  Continue
      If (iRow.GT.0) then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '*********** ERROR ***************'
         Write (Lu,*) 'CtoF and User-defined Coordinates'
         Write (Lu,*) 'are mutually exclusive.          '
         Write (Lu,*) '*********************************'
         Call Quit_OnUserError()
      EndIf
      iNull    = 0
      New_Line = 1
      lCtoF    = .True.
      Lu_UDIC  = 91
      FilNam='UDIC'
      call molcas_open(Lu_UDIC,FilNam)
      ReWind(Lu_UDIC)
      Line=Get_Ln(LuRd)
      Call UpCase(Line)
      Call FixEqualSign2(Line,LuRd,Lu_UDIC,iNull,New_Line)
      Write (Lu_UDIC,'(A)') Line
      Close(Lu_UDIC)
      Go To 999
*                                                                      *
****** CONS ************************************************************
*                                                                      *
*     Copy constraints definition into the UDC file, to be read
*     (after fixing and merging, if necessary) by DefInt2/Cllct2.
*
 9478 Continue
      If (.Not.Expert) Then
         Write (Lu,*)
         Write (Lu,*) ' ************ ERROR ***************'
         Write (Lu,*) ' Obsolete input standard!'
         Write (Lu,*) ' The CONSTRAINT section should'
         Write (Lu,*) ' be define in the &Gateway input.'
         Write (Lu,*)
         Write (Lu,*) ' To override add the EXPERT keyword'
         Write (Lu,*) ' to the top of the SLAPAF input.'
         Write (Lu,*) ' **********************************'
         Call Quit_OnUserError()
      End If
      If (External_UDC) Then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '****************** ERROR *********************'
         Write (Lu,*) 'Multiple definitions of constraints.'
         Write (Lu,*) 'Check Gateway and Slapaf inputs for conflicts!'
         Write (Lu,*) '**********************************************'
         Call Quit_OnUserError()
      End If
      Lu_UDC=20
      FilNam='UDC'
      Lu_UDC=IsFreeUnit(Lu_UDC)
      Call Molcas_Open(Lu_UDC,FilNam)
 318  Continue
      Line=Get_Ln(LuRd)
      Call UpCase(Line)
      Call LeftAd(Line)
      Write(Lu_UDC,'(A)') Trim(Line)
      If (Line(1:4).ne.'END') Go To 318
      Close(Lu_UDC)
      Go To 999
*                                                                      *
****** VDWB VdW correction both coordinate and Hessian *****************
*                                                                      *
981   iOptC = iOr(1024,iOptC)
      iOptC = iOr(2048,iOptC)
      Go To 999
*                                                                      *
****** NO VDWB VdW correction both coordinate and Hessian **************
*                                                                      *
984   Mask=iSetAll-2**10-2**11
      iOptC = iAnd(Mask,iOptC)
      Go To 999
*                                                                      *
****** VDWB VdW correction for coordinate only *************************
*                                                                      *
982   iOptC = iOr(2048,iOptC)
      Go To 999
*                                                                      *
****** NO VDWB VdW correction for coordinate only **********************
*                                                                      *
985   Mask=iSetAll-2**11
      iOptC = iAnd(Mask,iOptC)
      Go To 999
*                                                                      *
****** VDWB VdW correction for Hessian only ****************************
*                                                                      *
983   iOptC = iOr(1024,iOptC)
      Go To 999
*                                                                      *
****** NO VDWB VdW correction for Hessian only *************************
*                                                                      *
986   Mask=iSetAll-2**10
      iOptC = iAnd(Mask,iOptC)
      Go To 999
*                                                                      *
****** OLDF ************************************************************
*                                                                      *
903   lOld = .True.
      Go To 999
*                                                                      *
****** CART ************************************************************
*                                                                      *
918   CurviLinear = .False.
      Go To 999
*                                                                      *
****** THRS ************************************************************
*                                                                      *
*     read the gradient threshold
*
 908  Char=Get_Ln(LuRd)
      Call Get_F1(1,ThrEne)
      Call Get_F1(2,ThrGrd)
      ThrInp=.True.
      Go To 999
*                                                                      *
****** TOLE ************************************************************
*                                                                      *
*     read the constraints threshold
*
 909  Char=Get_Ln(LuRd)
      Call Get_F1(1,ThrCons)
      ThrCons=Abs(ThrCons)
      Go To 999
*                                                                      *
****** SUPS ************************************************************
*                                                                      *
*     Introduce supersymmetry
*     Input format
*     nsg                (number of super groups)
*     Reapeat nsg times
*     nmem, (ind.., i = 1, nmem)
*
 911  LSup = .True.
      Char=Get_Ln(LuRd)
      Call Get_I1(1,nSupSy)
      Call GetMem(' NSup ','Allo','Inte',ipNSup,NSUPSY)
      Call GetMem('iAtom ','Allo','Inte',ipAtom,nsAtom)
      iStrt = ipAtom
      Do 950 i = ipNSup, ipNSup+nSupSy-1
         Read(LuRd,*,Err=9630)iWork(i),
     &       (iWork(j),j=iStrt,iStrt+iWork(i)-1)
         iStrt = iStrt + iWork(i)
 950  Continue
      Go To 999
9630  Call WarningMessage(2,'Error in RdCtl_Slapaf')
      Write (Lu,*)
      Write (Lu,*) '************ ERROR ****************'
      Write (Lu,*) 'Error while reading supersymmetry.'
      Write (Lu,*) '***********************************'
      Call Quit_OnUserError()
*                                                                      *
****** HUPD ************************************************************
*                                                                      *
914   Char=Get_Ln(LuRd)
      Read(Char,*) Char
      Call UpCase(Char)
      If (Trim(Char).eq.'BFGS') Then
         iOptH = 4
c     Else If (Trim(Char).eq.'MEYER') Then
c        iOptH = iOr(1,iAnd(iOptH,32))
c     Else If (Trim(Char).eq.'BP') Then
c        iOptH = iOr(2,iAnd(iOptH,32))
      Else If (Trim(Char).eq.'NONE') Then
         iOptH = iOr(8,iAnd(iOptH,32))
      Else If (Trim(Char).eq.'MSP') Then
         iOptH = iOr(16,iAnd(iOptH,32))
      Else If (Trim(Char).eq.'EU') Then
         iOptH = iOr(64,iAnd(iOptH,32))
      Else If (Trim(Char).eq.'TS-BFGS') Then
         iOptH = iOr(128,iAnd(iOptH,32))
      Else
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '************ ERROR ****************'
         Write (Lu,*) 'Unsupported Hessian update method: ',Trim(Char)
         Write (Lu,*) '***********************************'
         Call Quit_OnUserError()
      End If
      Go To 999
*                                                                      *
****** MAXS ************************************************************
*                                                                      *
 915  Char=Get_Ln(LuRd)
      If (Char.eq.BLine) Go To 915
      If (Char(1:1).eq.'*') Go To 915
      Call Get_F1(1,Beta)
      Go To 999
*                                                                      *
****** PRIN ************************************************************
*                                                                      *
 920  Char=Get_Ln(LuRd)
      Call UpCase(Char)
      If (Char.eq.BLine) Go To 920
      If (Char(1:1).eq.'*') Go To 920
      Call Get_I1(1,mPrint)
      Do 921 i = 1, mPrint
 922     Char=Get_Ln(LuRd)
         Call UpCase(Char)
         If (Char.eq.BLine) Go To 922
         If (Char(1:1).eq.'*') Go To 922
         Call Get_I1(1,iRout)
         Call Get_I1(2,kPrint)
         nPrint(iRout)=kPrint
 921  Continue
      Go To 999
*                                                                      *
****** PRFC ************************************************************
*                                                                      *
*     set nPrint to print internal coordinates and hessian
*
9201  nPrint(21)=6  ! Eigen-/Values/Vectors of the Hessian (diagmtrx)
      nPrint(116)=6 ! Internal Forces (rlxctl)
      If (.NOT.Request_Alaska) nPrint(30)=6 ! Coord.s & Forces (defint)
      nPrint(122)=6 ! Auto-Defined Internal coordinates (printq_sl)
      Go To 999
*                                                                      *
****** ITER ************************************************************
*                                                                      *
*     read max iterations
*
 925  Char=Get_Ln(LuRd)
      Call Get_I1(1,iTmp)
      MxItr=Min(iTmp,MxItr)
      Go To 999
*                                                                      *
****** AI   ************************************************************
*                                                                      *
*     Activate Kriging
*
100   Char=Get_Ln(LuRd)
      If (Char.eq.'Kriging'.or.Char.eq.'kriging') then
       Kriging = .True.
      Else
       Call WarningMessage(1,'Illegal AI method selected.')
      EndIf
      Go To 999
*                                                                      *
****** AIAM ************************************************************
*                                                                      *
*      Analitical or numerical Mat'ern derivatives
*
101   Char=Get_Ln(LuRd)
      If (Char.eq.'False'.or.Char.eq.'false') anAI = .False.
      Go To 999
*                                                                      *
****** AIL  ************************************************************
*                                                                      *
*     Widht limits of the Mat`ern function
*
102   Char=Get_Ln(LuRd)
      Call Get_F(1,lb,3)
      Go To 999
*                                                                      *
****** AINX ************************************************************
*                                                                      *
*     The resolution of the predicted path
*
103   Char=Get_Ln(LuRd)
      Call Get_I(1,npxAI,1)
      Go To 999
*                                                                      *
****** AIP  ************************************************************
*                                                                      *
*     Parameter of differentiability for Mat`ern function
*
104   Char=Get_Ln(LuRd)
      Call Get_F(1,pAI,1)
      If(pAI.gt.2.or.pAI.lt.1) anAI = .False.
      Go To 999
*                                                                      *
****** AISP ************************************************************
*                                                                      *
*     Defining the number of source points for the AI method
*
105   Char=Get_Ln(LuRd)
      Call Get_I(1,nspAI,1)
      Go To 999
*                                                                      *
****** AIMI ************************************************************
*                                                                      *
*     Maximum number of Iterations for the AI method
*
106   Char=Get_Ln(LuRd)
      Call Get_I(1,miAI,1)
      Go To 999
*                                                                      *
****** AIME ************************************************************
*                                                                      *
*     Minimum egergy differences of the last two Iterations
*     (loop exit condition)
*
107   Char=Get_Ln(LuRd)
      Call Get_I(1,meAI,1)
      Go To 999
*                                                                      *
****** BAKE ************************************************************
*                                                                      *
926   Baker = .True.
      Go To 999
*                                                                      *
****** SCHL ************************************************************
*                                                                      *
927   Schlegel = .True.
      Go To 999
*                                                                      *
****** DDVS ************************************************************
*                                                                      *
9271  DDV_Schlegel = .True.
      Go To 999
*                                                                      *
****** NOLA ************************************************************
*                                                                      *
930   CallLast = .False.
      Go To 999
*                                                                      *
****** NOLI ************************************************************
*                                                                      *
928   Line_Search = .False.
      Go To 999
*                                                                      *
****** LAST ************************************************************
*                                                                      *
9280  Char=Get_Ln(LuRd)
      Call LeftAd(Char)
      If (Char.eq.BLine) Go To 9280
      If (Char(1:1).eq.'*') Go To 9280
      Call UpCase(Char)
      Call Put_cArray('LastEnergyMethod',Char,8)
      Go To 999
*                                                                      *
****** LINE ************************************************************
*                                                                      *
9281  Line_Search = .True.
      Go To 999
*                                                                      *
****** HWRS ************************************************************
*                                                                      *
929   HWRS=.True.
      Go To 999
*                                                                      *
****** WIND ************************************************************
*                                                                      *
 934  Char=Get_Ln(LuRd)
      Call UpCase(Char)
      If (Char.eq.BLine) Go To 934
      If (Char(1:1).eq.'*') Go To 934
      Call Get_I1(1,nWndw)
      Go To 999
*                                                                      *
****** NEWT ************************************************************
*                                                                      *
935   Mask = iSetAll
      Mask = Mask - 2**0 - 2**1 - 2**2 - 2**3
      iOptC = iOr(2**0,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** C1-D ************************************************************
*                                                                      *
936   Mask = iSetAll
      Mask = Mask - 2**0 - 2**1 - 2**2 - 2**3
      iOptC = iOr(2**1,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** C2-D ************************************************************
*                                                                      *
937   Mask = iSetAll
      Mask = Mask - 2**0 - 2**1 - 2**2 - 2**3
      iOptC = iOr(2**2,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** RATI ************************************************************
*                                                                      *
938   Mask = iSetAll
      Mask = Mask - 2**0 - 2**1 - 2**2 - 2**3
      iOptC = iOr(2**3,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** DXDX ************************************************************
*                                                                      *
939   Mask = iSetAll
      Mask = Mask - 2**4 - 2**5 - 2**6
      iOptC = iOr(2**4,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** DXG  ************************************************************
*                                                                      *
940   Mask = iSetAll
      Mask = Mask - 2**4 - 2**5 - 2**6
      iOptC = iOr(2**5,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** GG   ************************************************************
*                                                                      *
941   Mask = iSetAll
      Mask = Mask - 2**4 - 2**5 - 2**6
      iOptC = iOr(2**6,iAnd(iOptC,Mask))
      Go To 999
*                                                                      *
****** MODE ************************************************************
*                                                                      *
*-----Mode following algorithm
942   Continue
*     Read (5,'(A)',End=9610) Char
      Char=Get_Ln(LuRd)
      If (Char.eq.BLine) Go To 942
      If (Char(1:1).eq.'*') Go To 942
      Call Get_I1(1,mode)
      Go To 999
*                                                                      *
****** NUME ************************************************************
*                                                                      *
945   lNmHss = .True.
      Go To 999
*                                                                      *
****** THER ************************************************************
*                                                                      *
9451  lNmHss = .True.
      lTherm = .True.
      Char=Get_Ln(LuRd)
      Call Get_I1(1,nsRot)
      Char=Get_Ln(LuRd)
      Call Get_F1(1,UserP)
9454  Char=Get_Ln(LuRd)
      Call UpCase(Char)
      If (Char(1:4).eq.'END ') then
         If (nUserPT.EQ.0) then
           nUserPT=1
           UserT(1)=298.15d0
         EndIf
         Go To 999
      EndIf
      nUserPT=nUserPT+1
      Call Get_F1(1,UserT(nUserPT))
      Go To 9454
*                                                                      *
****** DISO ************************************************************
*                                                                      *
9452  lDoubleIso = .True.
      Go To 999
*                                                                      *
****** CUBI ************************************************************
*                                                                      *
947   Cubic  = .True.
      Go To 999
*                                                                      *
****** DELT ************************************************************
*                                                                      *
 946  Char=Get_Ln(LuRd)
      Call Get_F1(1,Delta)
      Go To 999
*                                                                      *
****** TS   ************************************************************
*                                                                      *
 951  Mask=iSetAll - 2**7
      iOptC=iAnd(Mask,iOptC)
      Go To 999
*                                                                      *
****** EXTR ************************************************************
*                                                                      *
*     Put the program name and the time stamp onto the extract file
*
971   Write (Lu,*)
     &'RdCtl_Slapaf: EXTRACT option is redundant and is ignored!'
      Go To 999
*                                                                      *
****** NOHW ************************************************************
*                                                                      *
 960  HWRS=.False.
      Go To 999
*                                                                      *
****** RTRN ************************************************************
*                                                                      *
 962  Char = Get_Ln(LuRd)
      Call UpCase(Char)
      Call Get_I1(1,Max_Center)
      Call Get_F1(2,rtrnc)
      If (Index(Char,'ANGSTROM').ne.0) Rtrnc = Rtrnc/angstr
      Go To 999
*                                                                      *
****** FIND ************************************************************
*                                                                      *
 963  FindTS=.True.
      Go To 999
*                                                                      *
****** TSCO ************************************************************
*                                                                      *
 320  LuTS=20
      FilNam='TSC'
      LuTS=IsFreeUnit(LuTS)
      Call Molcas_Open(LuTS,FilNam)
 319  Line=Get_Ln(LuRd)
      Call UpCase(Line)
      Call LeftAd(Line)
      Write(LuTS,'(A)') Trim(Line)
      If (Line(1:4).ne.'END') Go To 319
      Close(LuTS)
      TSConstraints=.True.
      Go To 999
*                                                                      *
****** FUZZ ************************************************************
*                                                                      *
 123  Char = Get_Ln(LuRd)
      Call UpCase(Char)
      Call Get_F1(1,rFuzz)
      If (Index(Char,'ANGSTROM').ne.0) rFuzz = rFuzz/angstr
      rFuzz=Max(rFuzz,1.0D-3)
      Go To 999
*                                                                      *
****** MEP-/MEP  *******************************************************
*                                                                      *
 964  MEP=.True.
      rMEP=.False.
      Go To 999
*                                                                      *
****** NMEP/NIRC *******************************************************
*                                                                      *
 965  Char=Get_Ln(LuRd)
      Call Get_I1(1,nMEP)
      nMEP=Min(Max(nMEP,1),MaxItr)
      Go To 999
*                                                                      *
****** MEPT/IRCT *******************************************************
*                                                                      *
 321  Char=Get_Ln(LuRd)
      Call UpCase(Char)
      If (Char(1:6).eq.'SPHERE') Then
         MEP_Type='SPHERE'
      Else If (Char(1:5).eq.'PLANE') Then
         MEP_Type='TRANSVERSE'
      Else
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '********** ERROR **********'
         Write (Lu,*) ' Unrecognized MEP/IRC type.'
         Write (Lu,*) '***************************'
         Call Quit_OnUserError()
      End If
      Go To 999
*                                                                      *
****** MEPA/IRCA *******************************************************
*                                                                      *
 322  Char=Get_Ln(LuRd)
      Call UpCase(Char)
      If (Char(1:2).eq.'GS') Then
         MEP_Algo='GS'
      Else If (Char(1:2).eq.'MB') Then
         MEP_Algo='MB'
      Else
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '************* ERROR ************'
         Write (Lu,*) ' Unrecognized MEP/IRC algorithm.'
         Write (Lu,*) '********************************'
         Call Quit_OnUserError()
      End If
      Go To 999
*                                                                      *
****** REFE ************************************************************
*                                                                      *
 966  Call GetMem('RefGeom','Allo','Real',ipRef,3*nsAtom)
      Call Read_v(LuRd,Work(ipRef),1,3*nsAtom,1,iErr)
      If (iErr.ne.0) Then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '************ ERROR ***************'
         Write (Lu,*) 'Error reading reference structure.'
         Write (Lu,*) '**********************************'
         Call Quit_OnUserError()
      End If
      Ref_Geom=.True.
      Go To 999
*                                                                      *
****** RS-P ************************************************************
*                                                                      *
 967  Mask=iSetAll
      Mask=Mask-2**9
      iOptC=iAnd(iOptC,Mask)
      Go To 999
*                                                                      *
****** GNRM ************************************************************
*                                                                      *
 968  Char = Get_Ln(LuRd)
      Call Get_F1(1,GNrm_Threshold)
      Go To 999
*                                                                      *
****** GRAD ************************************************************
*                                                                      *
 979  Call GetMem('ReGradient','Allo','Real',ipGradRef,nDimbc)

      Call Read_v(LuRd,Work(ipGradRef),1,nDimbc,1,iErr)
      If (iErr.ne.0) Then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '*************** ERROR ******************'
         Write (Lu,*) 'Error reading reference gradient vector.'
         Write (Lu,*) '****************************************'
         Call Quit_OnUserError()
      End If
*
*     If there is a transverse vector stored, we are not using this one
*
      Call qpg_dArray('Transverse',Found,nRP)
      If (.Not. Found) Ref_Grad=.True.
      Go To 999
*                                                                      *
****** rMEP ************************************************************
*                                                                      *
 980  rMEP=.True.
      MEP=.False.
      Go To 999
*                                                                      *
****** rHidden *********************************************************
*                                                                      *
 988  Line = Get_Ln(LuRd)
      Call UpCase(Line)
      Call Get_F1(1,rHidden)
      If (rHidden.lt.Zero) Then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '************ ERROR *****************'
         Write (Lu,*) 'Error reading rHidden. Should be >0.'
         Write (Lu,*) '************************************'
         Call Quit_OnUserError()
      End If
      If (Index(Line,'ANGSTROM').ne.0) rHidden = rHidden/angstr
      Go To 999
*                                                                      *
****** IRC *************************************************************
*                                                                      *
 997  Call Qpg_iScalar('IRC',Found)
      If (Found) Then
         Call Get_iScalar('IRC',IRC)
      Else
         IRC=1
         Call Put_iScalar('IRC',IRC)
      End If
      MEP=.True.
      rMEP=.False.
      Go To 999
*                                                                      *
****** MEPStep/IRCStep *************************************************
*                                                                      *
 9971 Char=Get_Ln(LuRd)
      Call UpCase(Char)
      Call Get_F1(1,dMEPStep)
*
*     Note that according to the Gonzalez-Schlegel method, only half
*     this step is used in the constraint
*
      If (Index(Char,'ANGSTROM').ne.0) dMEPStep = dMEPStep/angstr
      Go To 999
*                                                                      *
****** REAC ************************************************************
*                                                                      *
 996  Explicit_IRC=.True.
      Call GetMem('TmpRx','Allo','Real',ipTmpRx,3*nsAtom)
      Call Read_v(LuRd,Work(ipTmpRx),1,3*nsAtom,1,iErr)
      If (IErr.ne.0) Then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '********** ERROR ***********************'
         Write (Lu,*) ' Error while reading the Reaction vector'
         Write (Lu,*) '****************************************'
         Call Quit_OnUserError()
      End If
      Go To 999
*                                                                      *
****** HRMS ************************************************************
*                                                                      *
 995  HrmFrq_Show=.True.
      Go To 999
*                                                                      *
****** CNWE ************************************************************
*                                                                      *
 990  Char=Get_Ln(LuRd)
      Call Get_F1(1,CnstWght)
      Go To 999
*                                                                      *
****** NOEM ************************************************************
*                                                                      *
 991  eMEPTest=.False.
      Go To 999
*                                                                      *
****** EXPE ************************************************************
*                                                                      *
 993  Expert=.True.
      Go To 999
*                                                                      *
****** REDU ************************************************************
*                                                                      *
 994  Redundant=.True.
      Go To 999
*                                                                      *
****** FALC ************************************************************
*                                                                      *
 800  isFalcon=.True.
      Go To 999
*                                                                      *
****** TRAC ************************************************************
*                                                                      *
 910  Track=.True.
      Go To 999
*                                                                      *
************************************************************************
************************   End of input section   **********************
************************************************************************
*                                                                      *
 998  Continue
*                                                                      *
************************************************************************
*                                                                      *
*     Now start fixing constraints.
*     First get the external constraints.
*
      If (External_UDC) Then
        Call Merge_Constraints('UDC.Gateway','','UDC',nLambda,iRow_c)
      Else
        Call Merge_Constraints('','UDC','UDC',nLambda,iRow_c)
      End If
*
*     Initial preprocessing
*
      If (iRow_c.gt.1) Then
         Lu_UDC=IsFreeUnit(20)
         Call Molcas_Open(Lu_UDC,'UDC')
         Call Preprocess_UDC(Lu_UDC,iPrint)
         Close (Lu_UDC)
      Else
         NADC=.False.
      End If
*
*     Add NAC if needed
*
      If (NADC) Then
         Lu_UDCTMP=IsFreeUnit(20)
         Call Molcas_Open(Lu_UDCTMP,'UDCTMP')
         Write(Lu_UDCTMP,*) 'NADC = NAC'
         Write(Lu_UDCTMP,*) 'VALUE'
         Write(Lu_UDCTMP,*) 'NADC = 0.0'
         Write(Lu_UDCTMP,*) 'END'
         Close (Lu_UDCTMP)
         Call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
      End If
*
*     Add MEP/IRC if needed
*
      If (MEP.or.rMEP.or.(Abs(IRC).eq.1)) Then
         If (Abs(IRC).eq.1) Then
           MEPLab='IRC'
         Else
           MEPLab='MEP'
         End If
         If (MEPCons.and.(.Not.Expert)) Then
            Call WarningMessage(2,'Error in RdCtl_Slapaf')
            Write (Lu,*)
            Write (Lu,*)
     &'***************** ERROR ********************'
            Write (Lu,*)
     &' There is a '//Trim(Mep_Type)//' constraint that may'
            Write (Lu,*)
     &' conflict with '//MEPLab//' calculations.'
            Write (Lu,*)
     &' You should not explictly specify this constraint,'
            Write (Lu,*)
     &' but just rely on '//MEPLab//'Step/'//MEPLab//'Type keywords.'
            Write (Lu,*)
     &' If you really know what you are doing, you'
            Write (Lu,*)
     &' can use the EXPERT keyword.'
            Write (Lu,*)
     &'********************************************'
            Call Quit_OnUserError()
         End If
         WeightedConstraints=.True.
         Valu=dMEPStep
         If (MEP_Type.eq.'SPHERE') Valu=Abs(Valu)
         If (MEP.and.(MEP_Algo.eq.'GS')) Valu=Half*Valu
         If (rMEP) Valu=Max(Dble(iMEP+1),One)*Valu
         Lu_UDCTMP=IsFreeUnit(20)
         Call Molcas_Open(Lu_UDCTMP,'UDCTMP')
         Write(Lu_UDCTMP,*) MEPLab//' = '//MEP_Type
         Write(Lu_UDCTMP,*) 'VALUE'
         Write(Lu_UDCTMP,*) MEPLab//' = ',Valu
         Write(Lu_UDCTMP,*) 'END'
         Close (Lu_UDCTMP)
         Call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
         Beta=Abs(dMEPStep)
         MEPnum=nLambda
      End If
*
*     Final fixes
*
      Call Fix_UDC(iRow_c,nLambda,AtomLbl,nsAtom,nStab,.True.)
*                                                                      *
************************************************************************
*                                                                      *
*     Initiate some variables which can only be set after the input has
*     been read.
*
      If ((.Not.ThrInp).and.(.Not.Baker)) ThrEne=Zero
      Call Init2
*     Gradients are not needed at the first iteration of a numerical
*     Hessian procedure (and only that, i.e. MxItr=0)
      FirstNum = (lRowH.or.lNmHss.or.Cubic)
     &           .and.(Iter.eq.1).and.(MxItr.eq.0)
      If ((SuperName.eq.'slapaf').and.(.not.FirstNum)) Then
         If (Track) Then
            Call Process_Track()
         Else
            Call Put_iArray('Root Mapping',RootMap,0)
         End If
         Call Process_Gradients()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put in the "Reaction vector" in Cartesians.
*     Priority order:
*     1) Explicit by user input (REAC keyword)
*     2) Found on RunOld
*     3) Found on RunFile
*
      If (Abs(IRC).eq.1) Then
*
*        If this is the first macro iteration in the IRC search then
*        pick up the reaction vector.
*
         If (Explicit_IRC.and.iMEP.eq.0) Then
*           Case 1)
            call dcopy_(3*nsAtom,Work(ipTmpRx),1,Work(ipMF),1)
         Else If (iMEP.eq.0) Then
            Call NameRun('RUNOLD')
            Call qpg_dArray('Reaction Vector',Found,nRx)
C           Write (6,*) 'RUNOLD: Found=',Found
            If (Found) Then
*              Case 2)
               Call Get_dArray('Reaction Vector',Work(ipMF),3*nsAtom)
               Call NameRun('RUNFILE')
            Else
               Call NameRun('RUNFILE')
               Call qpg_dArray('Reaction Vector',Found,nRx)
C              Write (6,*) 'RUNFILE: Found=',Found
               If (Found) Then
*                 Case 3)
                  Call Get_dArray('Reaction Vector',Work(ipMF),3*nsAtom)
               Else
                  Call WarningMessage(2,'Error in RdCtl_Slapaf')
                  Write (6,*)
                  Write (6,*) '************ ERROR **************'
                  Write (6,*) 'IRC calculation but no IRC vector'
                  Write (6,*) '*********************************'
                  Call Quit_OnUserError()
               End If
            End If
         End If
*
*        Fix the direction forward/backwards
*
         If (iMEP.eq.0.and.iRC.eq.-1) Call DScal_(3*nsAtom,-1.0D0,
     &                                           Work(ipMF),1)
         If (iMEP.eq.0.and.MEP_Type.eq.'TRANSVERSE')
     &      Call Put_dArray('Transverse',Work(ipMF),3*nsAtom)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (FindTS.and.(.Not.TSConstraints)) Then
         Call SysWarnMsg('RdCtl_Slapaf','WARNING:',
     &   'FindTS specified, but no TSConstraints. '//
     &   'It is highly recommended to use TSConstraints in SLAPAF '//
     &   'instead of (or in addition to) global constraints when '//
     &   'using FindTS. '//
     &   'TSConstraints will be lifted in the final TS search.')
      End If
      TSConstraints=TSConstraints.and.FindTS
*
      If ((MEP.or.rMEP).and.(.NOT.Request_Alaska)) Then
*
*        If no initial direction given, use the gradient (force)
*
         Call qpg_dArray('Transverse',Found,nRP)
         If (.Not.Found.And..Not.Ref_Grad) Then
*        Assume the initial reaction vector is already massaged
            If (Explicit_IRC) Then
               Call Put_dArray('Transverse',Work(ipTmpRx),3*nsAtom)
            Else
*        The direction is given by the gradient, but in weighted coordinates
               Call Allocate_Work(ipDir,3*nsAtom)
               iOff=0
               ip_Grd=ipGx+(iter-1)*3*nsAtom
               Do iAtom=1,nsAtom
                  xWeight=Work(ipWeights+iAtom-1)
                  Do ixyz=1,3
                     Work(ipDir+iOff)=Work(ip_Grd+iOff)/xWeight
                     iOff=iOff+1
                  End Do
               End Do
               Call Put_dArray('Transverse',Work(ipDir),3*nsAtom)
               Call Free_Work(ipDir)
            End If
         End If
      End If
*
      If (Explicit_IRC) Call Free_Work(ipTmpRx)
*
*     Activate MPS update of Hessian if FindTS
*
      If (FindTS) Then
*
         If (iAnd(iOptH,64).eq.64) Then
            iOptH=iOr(64,iAnd(iOptH,32)) ! EU
         Else If (iAnd(iOptH,128).eq.128) Then
            iOptH=iOr(128,iAnd(iOptH,32)) ! TS-BFGS
         Else
            iOptH=iOr(16,iAnd(iOptH,32)) ! MSP
         End If
         iOptC=iOr(iOptC,4096)
*
*------- Increase the update window so that we will not lose the update
*        which generated the negative curvature.
*
         nWndw=4*nWndw
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Modify some options if TS search
*
      If (iAnd(iOptC,128).ne.128) Then
         If (iAnd(iOptH,8).ne.8) iOptH=iOr(16,iAnd(iOptH,32)) ! MSP
         ThrB=0.01D+00
         Line_search=.False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     For TS optimization with the Saddle method set update to MSP
*
      Call qpg_dArray('Saddle',Found,nSaddle)
      If (Found.and.nSaddle.ne.0) Then
         Call Allocate_Work(ipTmp,nSaddle)
         Call Get_dArray('Saddle',Work(ipTmp),nSaddle)
         HSR0=Work(ipTmp+nSaddle-3)
         HSR=Work(ipTmp+nSaddle-2)
         Update=Work(ipTmp+nSaddle-1)
         If (Update.eq.2.0d0) Then
*
*           Enable FindTS procedure
*
C           Write (6,*) 'Enable FindTS procedure'
            If (iAnd(iOptH,8).ne.8) iOptH=iOr(16,iAnd(iOptH,32)) ! MSP
            nWndw=4*nWndw
*           make it look as if this were FindTS with constraints
            FindTS=.True.
            TSConstraints=.True.
            iOptC=iOr(iOptC,4096)
            iOptC=iOr(iOptC,8192)
            Beta=0.1d0
*
         Else
*
*           Normal constrained optimization with a reduced threshold.
*           Let the threshold be somewhat tighter as we are close to
*           the TS.
*
            If (HSR/HSR0.lt.0.20D0.or.HSR.lt.0.20D0) Then
*              ThrGrd=0.0003D0
               Beta=0.1d0
            Else
*              ThrGrd=0.003D0
               ThrGrd=Ten*ThrGrd
            End If
*
*           Add the constraints from the Saddle method
*
            Call Merge_Constraints('UDC','UDC.Saddle','UDC',
     &                             nLambda,iRow_c)

         End If
         Call Free_Work(ipTmp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Modify some options if constraints are part of the calculation.
*
      If ((nLambda.gt.0).or.TSConstraints) Then
         iOptC=iOr(iOptC,256) ! Constraints
         Line_search=.False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     No iterations set iOptC=0
*
      If (MxItr.eq.0) iOptC=0
*                                                                      *
************************************************************************
*                                                                      *
*     Activate some additional printing for numerical Hessian
*
CGGd: Coherency with patch 7.1.615 !      If (lNmHss) nPrint(122)=10
*                                                                      *
************************************************************************
*                                                                      *
*.....Do some preprocessing due to input choice
*
      If (Request_Alaska) nPrint(51)=0
      Call PrePro(iRow,iInt,nFix,nsAtom,mInt,Work(ipCoor))
*                                                                      *
************************************************************************
*                                                                      *
*.....Write out input parameters, No output if we didn't find the proper
*     gradients on the runfile. We will come back!!!
*
      If (SuperName.eq.'slapaf') Then
         If (.Not.Request_Alaska) Call WrInp_sl(iRow)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      User_Def = iRow.ne.0
      If (.Not.User_Def) nBVec=1
*                                                                      *
************************************************************************
*                                                                      *
      If (lNmHss.and.lRowH) then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
       Write (Lu,*)
       Write (Lu,*) '**************************************************'
       Write (Lu,*) ' ERROR: NUMErical and ROWH are mutually exclusive '
       Write (Lu,*) '**************************************************'
       Call Quit_OnUserError()
      EndIf
      If (lCtoF.and.User_Def) then
         Call WarningMessage(2,'Error in RdCtl_Slapaf')
         Write (Lu,*)
         Write (Lu,*) '******************************************'
         Write (Lu,*) ' ERROR: CtoF and User-defined Coordinates '
         Write (Lu,*) '        are mutually exclusive.           '
         Write (Lu,*) '******************************************'
         Call Quit_OnUserError()
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('RdCtl_Slapaf')
      Return
*
      End
      Subroutine FixEqualSign(Line,LuRd)
      Character*(*) Line
      Character*180 Temp_Line
      Character*180 Get_Ln
      External Get_Ln
*
      nLine=LEN(Line)
      If (nLine.gt.LEN(Temp_Line)) Then
         Call WarningMessage(2,'Error in FixEqualSign!')
         Call Abend()
      End If
*
      Temp_Line=Line
      Call LeftAd(Temp_Line)
      ix = iCLast(Temp_Line,nLine)
      Temp_Line(ix+2:ix+2)='='
      ix = ix + 2
*
      Line=Get_Ln(LuRd)
      Call LeftAd(Line)
      iy = iCLast(Line,nLine)
      If (ix+2+iy.gt.nLine) Then
         Call WarningMessage(2,'Problems merging lines!')
         Call Abend()
      End If
      Temp_Line(ix+2:nLine) = Line(1:nLine-ix-1)
      Line = Temp_Line
      Call UpCase(Line)
*
      Return
      End
*
      Subroutine FixEqualSign2(Line,LuRd,Lu_UDIC,iRow,New_Line)
      Character*(*) Line
      Character*180 Temp_Line
      Character*180 Get_Ln
      External Get_Ln
*
      nLine=LEN(Line)
      If (nLine.gt.LEN(Temp_Line)) Then
         Call WarningMessage(2,'Error in FixEqualSign!')
         Call Abend()
      End If
*
      Temp_Line=Line
      Call LeftAd(Temp_Line)
      ix = iCLast(Temp_Line,nLine)
*                                                                      *
************************************************************************
*                                                                      *
*     Read the next line and determine if the lines should be merged.
*
      Line=Get_Ln(LuRd)
      Call LeftAd(Line)
      iy = iCLast(Line,nLine)
      Call UpCase(Line)
*                                                                      *
************************************************************************
*                                                                      *
      If (Index(Line(1:iy),'END ').eq.1) Then
         iRow = iRow + 1
         Write (Lu_UDIC,'(A)') Temp_Line
         New_Line=2
*                                                                      *
************************************************************************
*                                                                      *
*     If the line contains two or more items we should merge the lines.
*
      Else If (Index(Line(1:iy),' ').eq.0) Then
*
*        Just one item
*
         iRow = iRow + 1
         Write (Lu_UDIC,'(A)') Temp_Line
         New_Line=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Temp_Line(ix+2:ix+2)='='
         ix = ix + 2
         If (ix+2+iy.gt.nLine) Then
            Call WarningMessage(2,'Problems merging lines!')
            Call Abend()
         End If
         Temp_Line(ix+2:nLine) = Line(1:nLine-ix-1)
         Line = Temp_Line
         Call UpCase(Line)
         New_Line=0
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
