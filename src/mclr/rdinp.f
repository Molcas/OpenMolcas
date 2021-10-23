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
       Subroutine RdInp_MCLR
************************************************************************
*                                                                      *
*     Locate input stream and read commands                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"
#include "negpre.fh"
#include "sa.fh"
      Parameter ( nCom=38 )
      Character*72 Line
      Character*4 Command,ComTab(nCom)
      Character*8 Label
      Character*2 Element(MxAtom)
      Logical     Epsilon_Undef
#include "chomclr.fh"
      Data ComTab/'TITL','DEBU','ROOT','EXTR','PRCI',
     &            'PROR','ITER','THRE','END ','TIME',
     &            'CALC','NOFI','SEWA','NOCO','NOTW',
     %            'SPIN','PRIN','PCGD','RESI','NOTO',
     &            'EXPD','NEGP','LOWM','ELHE','SAVE',
     &            'RASS','DISO','CASI','SALA','NODE',
     &            'ESTE','MOUT','MASS','NAC ','$$$$',
     &            'THER','NEWC','TWOS'/
      Dimension idum(1)
*----------------------------------------------------------------------*
*     Locate "start of input"                                          *
*----------------------------------------------------------------------*
      Call RdNLst(5,'MCLR')
*----------------------------------------------------------------------*
*     Define default values                                            *
*----------------------------------------------------------------------*
      debug=.False.
      Epsilon_Undef=.True.
      Call Get_info_Static(iDum(1))
      istate=1     ! State for which the Lagrangian is calc.
      override=.false.
      If (debug) write(6,*) 'Got info.fh'
      lRoots=-1
      kprint=0
      ngp=.false.
      NoFile=.false.
      mTit=0
      Omega=0.0d0
      elechess=.false.
      TimeDep=.false.
      PrCI=.false.
      nexp_max=100
      CIthrs=0.05d0
      PrOrb=.false.
      SewLab='NONE    '
      Page=.false.
      OEthrs=1.0d0
      ibreak=2
      nIter=200
      RASSI=.false.
      spinpol=.false.
      SA=.false.
      esterr=.false.
      FANCY_PRECONDITIONER=.true.
      newpre=.true.
      save=.false.
      isotop=.true.
      Call lCopy(mxAtm*3+3,[.true.],0,lCalc,1)
      Do i=1,nDisp
       DspVec(i)=i
      End Do
      nmode = 0
      lmass = .false.
      CASINT=.true.
      NACstates(1)=0
      NACstates(2)=0
      NSSA(1)=0
      NSSA(2)=0
      isNAC=.false.
      NewCho=.false.
*Cholesky. Cannot modify it in the input (yet?)
      dmpk=1.0d-2
      Nscreen=10
      Deco=.true.
      Update=.true.
      Estimate=.false.
      TwoStep=.false.
      StepType='xxxx'
*----------------------------------------------------------------------*
*     Read the input stream line by line and identify key command      *
*----------------------------------------------------------------------*
100   Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 100
      Call StdFmt(Line,Command)
      jCom=0
      Do iCom=1,nCom
        If ( Command.eq.ComTab(iCom) ) jCom=iCom
      End Do
      If ( jCom.eq.0 ) Then
         Write (6,'(A,A)') 'RdInp: illegal command:',Command
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Branch to the processing of the command sections                 *
*----------------------------------------------------------------------*
110   Select Case (jCom)
        Case (1)
          Go to 10
        Case (2)
          Go to 16
        Case (3)
          Go to 20
        Case (4)
          Go to 30
        Case (5)
          Go to 40
        Case (6)
          Go to 50
        Case (7)
          Go to 60
        Case (8)
          Go to 70
        Case (9)
          Go to 99
        Case (10)
          Go to 80
        Case (11)
          Go to 55
        Case (12)
          Go to 175
        Case (13)
          Go to 185
        Case (14)
          Go to 165
        Case (15)
          Go to 166
        Case (16)
          Go to 177
        Case (17)
          Go to 178
        Case (18)
          Go to 179
        Case (19)
          Go to 180
        Case (20)
          Go to 191
        Case (21)
          Go to 192
        Case (22)
          Go to 193
        Case (23)
          Go to 194
        Case (24)
          Go to 195
        Case (25)
          Go to 196
        Case (26)
          Go to 788
        Case (27)
          Go to 789
        Case (28)
          Go to 198
        Case (29)
          Go to 199
        Case (30)
          Go to 200
        Case (31)
          Go to 201
        Case (32)
          Go to 202
        Case (33)
          Go to 203
        Case (34)
          Go to 204
        Case (35)
          Go to 205
        Case (36)
          Go to 206
        Case (37)
          Go to 210
        Case (38)
          Go to 220
      End Select
*---  TITL ------------------------------------------------------------*
10    Continue
15    Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.'*' ) Goto 15
      Call StdFmt(Line,Command)
      jCom=0
      Do iCom=1,nCom
        If ( Command.eq.ComTab(iCom) ) jCom=iCom
      End Do
      If ( jCom.ne.0 ) Goto 110
      mTit=mTit+1
      If ( mTit.le.mxTit ) then
         Read(Line,'(18A4)') (TitleIN(iTit),iTit=(mTit-1)*18+1,mTit*18)
         Goto 15
      End If
      Goto 100
*
*----      ------------------------------------------------------------
788   RASSI=.true.
      If (debug) Write(6,*) 'Output for RASSI'
      goto 100
*----      ------------------------------------------------------------
789   double=.true.    ! Make double isotope substitutions
      goto 100
*---- DEBU ------------------------------------------------------------
16    debug=.true.
      Goto 100
*----      ------------------------------------------------------------
195   elechess=.true.
      If (debug) Write(6,*) 'Electric response'
      goto 100
*---- LOWM ------------------------------------------------------------
194   page=.true.
      If (debug) Write(6,*) 'Page memory'
      goto 100
*----      ------------------------------------------------------------
191   newpre=.false.
      If (debug) Write(6,*) 'New conditioner'
      goto 100
*----      ------------------------------------------------------------
196   SAVE=.TRUE.
      If (debug) Write(6,*) 'old integrals, not supported'
      goto 100
*----      ------------------------------------------------------------
198   CASINT=.true.
      If (debug) Write(6,*) 'CASPT2 integrals'
      goto 100
*---- EXPD ------------------------------------------------------------
192   Read(5,*) nexp_max
      If (debug) Write(6,*) 'Maximum explicit preconditioner',
     &           nexp_max
      goto 100
*----      ------------------------------------------------------------
179   iBreak=1
      Read(5,*) epsilon
      Epsilon_Undef=.False.
      If (debug) Write(6,*) 'Threshold:',epsilon
      Goto 100
*----      ------------------------------------------------------------
180   iBreak=2
      Read(5,*) epsilon
      Epsilon_Undef=.False.
      If (debug) Write(6,*) 'Threshold:',epsilon
      Goto 100
*----      ------------------------------------------------------------
178   Read(5,*) kprint
      If (debug) Write(6,*) 'Print level: ',kprint
      goto 100
*----      ------------------------------------------------------------
193   NGP=.true.
      If (debug) Write (6,*) 'NGP set to true'
      goto 100
*---- SALA ------------------------------------------------------------
199   SA=.true.
      Read(5,*) istate
      override=.true.
      If (debug) Write(6,*) 'Lagrangian for state: ',istate
      goto 100
*----      ------------------------------------------------------------
201   esterr=.true.
      goto 100
*---- NODE ------------------------------------------------------------
200   FANCY_PRECONDITIONER=.false.
      If (debug) Write(6,*) 'Turned of the fancy pcg'
      goto 100
*----      ------------------------------------------------------------
177   SPINPOL=.true.
      ispop=1
      If (debug) Write(6,*) 'RHF lagrangian, not supported'
      goto 100
*----      ------------------------------------------------------------
166   Do i=1,nDisp
       NTPert(i)=iAnd(nTPert(i),247)
      End Do
      Goto 100
*----      ------------------------------------------------------------
165   Do i=1,nDisp
       NTPert(i)=iAnd(nTPert(i),251)
      End Do
      Goto 100
*----      ------------------------------------------------------------

 185  Continue
 186  Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      Call StdFmt(Line,Command)
      If ( Command(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 186
      If (debug) Write(6,*) 'SEWARD INPUT'
      If ( Command(1:4).eq.'END '.or.Command(1:4).eq.'ENDS' ) Goto 100
      Read(Line,'(A8,I2,I2)',Err=998,End=999) SewLab,isym,ip
      iRc=-1
      iOpt=1
      iComp=ip
      iSyLbl=2**isym
      Label=SewLab
      Call iRdOne(iRc,iOpt,Label,iComp,idum,iSyLbl)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'RdInp: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call QTrace
         Call Abend()
      End If
*
*---  read number of symm. species ------------------------------------*
*
      ipp=0
      Do is=isym+1,nsym
       ipp=ipp+ldisp(is)
      end do
      DO id=nDisp,ndisp-ipp+1,-1
       DspVec(id+1)=dspVec(id)
       ntpert(id+1)=ntpert(id)
       lcalc(id+1)=lcalc(id)
      End Do
      id=ndisp-ipp+1
      DspVec(id)=ip
      ldisp(isym)=ldisp(isym)+1
      ndisp=ndisp+1
      ntpert(id)=2
      lcalc(id)=.true.
      SwLbl(id)=SewLab
      goto 185

 175  Nofile=.true.
      If (debug) Write(6,*) 'NOFILE      '
      Goto 100

*---  Process the "root" input card -----------------------------------*
20    Continue
25    Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 25
      Read(Line,*,Err=998,End=999) lRoots
      If (debug) Write(6,*) 'LROOT'
      Goto 100
*
*---  CALC ------------------------------------------------------------*
55    Continue
      Write (6,*) 'CALC is disabled!'
      Goto 100
*---  Process the "extract" input card --------------------------------*
30    Write (6,*) 'RdInp: EXTRACT option is redundant and is ignored!'
      Goto 100
*---  Process the "PrCI" input card -----------------------------------*
40    Continue
      PrCI=.true.
45    Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 45
      Read(Line,*,Err=998,End=999) CIthrs
      Goto 100
*---  Process the "PrOr" input card -----------------------------------*
50    Continue
      PrOrb=.true.
      Goto 100
*---  Process the "ITER" input card -----------------------------------*
60     Continue
65    Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 65
      Read(Line,*,Err=998,End=999) nIter
      Goto 100
*---  Process the "THRE" input card -----------------------------------*
70    Continue
75    Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 75
      Read(Line,*,Err=998,End=999) epsilon
      Epsilon_Undef=.False.
      Goto 100
*---  Process the "TIME" input card -----------------------------------*
80    Continue
      Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 80
      Read(Line,*,Err=998,End=999) Omega
      TimeDep=.true.
      nIter=100
      Goto 100
*---  Process the "MOUT" input card -----------------------------------*
202   Continue
      Write (6,*) 'MOUT is disabled!'
      Goto 100
*---  Process the "MASS" input card -----------------------------------*
203   Continue
      lmass = .true.
      iMass = 0
      Call Get_Name_All(Element)
*
*     Find out how many different elements are present in the molecule.
*
      Do i = 1, nAtoms
         If (Element(i).ne.'  ') iMass = iMass + 1
         Do j = i+1, nAtoms
            If (Element(j).eq.Element(i)) Element(j)='  '
         End Do
      End Do
      Do i=1,iMass
         Read(5,'(A3)')    cmass(i)
         Read(5,'(F15.8)') umass(i)
      End Do
*
*     Put the Info on the run file.

      Call Put_iScalar('iMass',iMass)
      Call Put_cArray('cmass',cmass(1),3*iMass)
      Call Put_dArray('umass',umass,iMass)
*
      Goto 100
*---  Process the "NAC " input card -----------------------------------*
204   Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' .or. Line(1:1).eq.'*' ) Goto 25
      Read(Line,*,Err=998,End=999) NACstates(1),NACstates(2)
      isNAC=.true.
      override=.true.
      If (debug) Write(6,*) 'Non-adiabatic couplings for states: ',
     &           NACstates(1),NACstates(2)
      Goto 100
*---  Process the "$$$$" input card -----------------------------------*
205   continue
*     not used
      Goto 100
*---  Process the "THERmochemistry input card -------------------------*
206   Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If (Line(1:1).eq.'*' ) Goto 206
      Read(Line,*,Err=998,End=999) nsRot
2060  Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If (Line(1:1).eq.'*' ) Goto 2060
      Read(Line,*,Err=998,End=999) UserP
2061  Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      If (Line(1:1).eq.'*' ) Goto 2061
      Call UpCase(Line)
      If (Line(1:4).eq.'END ') then
         If (nUserPT.EQ.0) then
           nUserPT=1
           UserT(1)=298.15d0
         EndIf
         GoTo 100
      EndIf
      nUserPT=nUserPT+1
      Read(Line,*,Err=998,End=999) UserT(nUserPT)
      Goto 2061
*---  Process the "NEWCho input card ----------------------------------*
210   NewCho=.True.
      Goto 100
*---  Process the "TWOStep" input card --------------------------------*
220   Read(5,'(A)',Err=998,End=999) Line
      Call LeftAd(Line)
      Call UpCase(Line)
      If (Line(1:1).eq.'*' ) Goto 220
      Read(Line,*,Err=998,End=999) StepType
      If (debug) Write(6,*) 'TWOSTEP kind: '//StepType
      If((StepType(1:4).ne.'FIRS').and.(StepType(1:4).ne.'SECO').and.
     &   (StepType(1:4).ne.'RUN1').and.(StepType(1:4).ne.'RUN2')) Then
         Call WarningMessage(2,'TWOStep: input error!')
         Call Quit_OnUserError()
      End If
      If (StepType(1:4).eq.'FIRS') then
        StepType(1:4)='RUN1'
      End If
      If (StepType(1:4).eq.'SECO') then
        StepType(1:4)='RUN2'
      End If
      TwoStep=.true.
      If (debug) Write(6,*) 'TWOSTEP kind: '//StepType
      Goto 100
*----------------------------------------------------------------------*
*     "End of input"                                                   *
*----------------------------------------------------------------------*
99    Continue
      Do i=1,3
         isym=irrfnc(2**(i-1))+1
         ipp=0
         Do is=isym+1,nsym
            ipp=ipp+ldisp(is)
         End Do
         Do id=nDisp,ndisp-ipp+1,-1
            DspVec(id+1)=dspVec(id)
            ntpert(id+1)=ntpert(id)
*           lcalc(id+1)=lcalc(id)
            Swlbl(id+1)=Swlbl(id)
         End Do
         id=ndisp-ipp+1
         DspVec(id)=i
         ldisp(isym)=ldisp(isym)+1
         ndisp=ndisp+1
         ntpert(id)=2
*        lcalc(id)=.true.
         write (Swlbl(id),'(a,i2)') 'MLTPL ',1
         iRc=-1
         iOpt=1
         Call iRdOne(iRc,iOpt,swlbl(id),dspvec(id),idum,iSyLbl)
      End Do
*
      If (Timedep) Then
         Do i=1,ndisp
            ntpert(i) = ior(ntpert(i),32)
         End Do
      End if
*
      If (Epsilon_Undef) Then
*        If (SA) Then
*           Epsilon=1.0D-6
*        Else
            Epsilon=1.0D-4
!        This I need to change back
*        End If
      End If
*
      If (debug) Write(6,*) 'FINITO'
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
*
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
 998  Write (6,*) 'RdInp: Error while reading input'
      Write (6,'(A,A)') 'Last command:',Line
      Call QTrace
      Call Abend()
999   Write (6,*) 'RdInp: Premature end of input file'
      Write (6,'(A,A)') 'Last command:',Line
      Call QTrace
      Call Abend()
      End
