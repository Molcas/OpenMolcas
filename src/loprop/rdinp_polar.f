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
      SubRoutine RdInp_Polar(LuSpool,NoField,Delta,MpProp_Level,
     &                       Bond_Threshold,iPlot,iPrint,Standard,
     &                       Opt_Method,UserDen,PrintDen,SubtractDen,
     &                       SubScale,Restart,TDensity,nStateI,nStateF,
     &                       XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,
     &                       ThrsMul,Alpha,LIonize)
*
      Implicit Real*8 (a-h,o-z)
*
*---- Define local variables
      Character*180 Key, Line
      Character*180 Get_Ln
      Character*12  Opt_Method
      External Get_Ln
      Logical NoField, Standard, UserDen, PrintDen, SubtractDen
      Logical Restart, Found, TDensity, XHole, Diffuse(3)
      Logical LIonize
      Dimension dLimmo(2)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Delta=0.001D00
      Call Get_iScalar('Highest Mltpl',lMax)
      MpProp_Level = lMax
      Bond_Threshold = 1.5D0
      iPlot = 0
      iPrint = 0
      Standard= .True.
      UserDen=.False.
      PrintDen=.False.
      SubtractDen=.False.
      Restart=.False.
      iRestart = 0
      SubScale=1.0D0
      nStateI=1
      nStateF=2
      Opt_Method = ' '
      TDensity=.False.
      XHole=.False.
      Diffuse(1)=.False.
      Diffuse(2)=.False.
      Diffuse(3)=.False.
      dLimmo(1)=0.65d0
      dLimmo(2)=2.0d0
      Thrs1=1d-5
      Thrs2=1d-4
      nThrs=3
      ThrsMul=1d-2
      Alpha=7.1421297D0
      LIonize=.False.

*     Comment on Alpha:
*     This value of Alpha is for backward-compability.
*     For large systems it may have to be reduced, for
*     example to 2.0.

*
*---- Locate "start of input"
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'LoProp')
*
  999 Continue
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)
*
c     If (Line(1:4).eq.'TITL') Go To 8000
      If (Line(1:4).eq.'NOFI') Go To 8001
      If (Line(1:4).eq.'DELT') Go To 8002
      If (Line(1:4).eq.'EXPA') Go To 8003
      If (Line(1:4).eq.'MPPR') Go To 8004
      If (Line(1:4).eq.'BOND') Go To 8005
      If (Line(1:4).eq.'PLOT') Go To 8006
      If (Line(1:4).eq.'PRIN') Go To 8007
      If (Line(1:4).eq.'USER') Go To 8008
      If (Line(1:4).eq.'PRDE') Go To 8009
      If (Line(1:4).eq.'SUBD') Go To 8010
      If (Line(1:4).eq.'REST') Go To 8011
      If (Line(1:4).eq.'TDEN') Go To 8012
      If (Line(1:4).eq.'XHOL') Go To 8013
      If (Line(1:4).eq.'DIFF') Go To 8014
      If (Line(1:4).eq.'ALPH') Go To 8015
      If (Line(1:4).eq.'LION') Go To 8016
      If (Line(1:4).eq.'END ') Go To 9000
      Write (6,*) 'Unidentified key word:', Key
      Call FindErrorLine
      Call Quit_OnUserError()
*
*>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C8000 Continue
      GoTo 999
*
*>>>>>>>>>>>>> NOFI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8001 Continue
      NoField=.True.
      GoTo 999
*
*>>>>>>>>>>>>> DELT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8002 Continue
      Key = Get_Ln(LuSpool)
      Call Get_F1(1,Delta)
      GoTo 999
*
*>>>>>>>>>>>>> EXPA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8003 Continue
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)

      If (Line(1:4) .eq. 'MIDP') Then
         Standard=.True.
      Else If (Line(1:4) .eq. 'OPTI') Then
         Standard=.False.
         Opt_Method = 'Optimized'
      Else If (Line(1:4) .eq. 'MULT') Then
         Standard = .False.
         Opt_Method = 'Multipole'
      Else
         Write(6,*) 'Undefined option for ''EXPAnsion center'':',Key
         Call FindErrorLine
         Call Quit_OnUserError()
      End If
      GoTo 999
*
*>>>>>>>>>>>>> MPPR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Read max multipole level for output in the MpProp file
 8004 Continue
      Key = Get_Ln(LuSpool)
      Call Get_I1(1,MpProp_Level)
      GoTo 999
*
*>>>>>>>>>>>>> BOND <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Read max bond length - all bonds longer than this will be ignored
 8005 Continue
      Key = Get_Ln(LuSpool)
      Call Get_F1(1,Bond_Threshold)
      GoTo 999
*
*>>>>>>>>>>>>> PLOT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Provide informations needed for plotting t vs. bond coordinate
* as well as printing a table of t values.
 8006 Continue
      iPlot = 1
      Goto 999
*
*>>>>>>>>>>>>> PRIN   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Read print level
 8007 Continue
      Key = Get_Ln(LuSpool)
      Call Get_I1(1,iPrint)
      GoTo 999
*
*>>>>>>>>>>>>>> USER   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Tell LoProp to read in density matrix(es) supplied by user.
 8008 Continue
      UserDen=.True.
      GoTo 999
*
*
*>>>>>>>>>>>>>> PRDE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Tell LoProp to output density matrix to file
 8009 Continue
      PrintDen=.True.
      NoField=.True. !Only static properties allowed
      GoTo 999
*
*
*>>>>>>>>>>>>>> SUBD   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Tell LoProp to subtract a user supplied density matrix
* from the current one (which could also be a USER supplied
* density matrix). Also input a scaling factor which the
* difference density is multiplied by.
 8010 Continue
      SubtractDen=.True.
      NoField=.True. !Only static properties allowed
      Key = Get_Ln(LuSpool)
      Call Get_F1(1,SubScale)
      GoTo 999
*
*>>>>>>>>>>>>> NOFI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Restart LoProp with the info stored in a previous LoProp
* calculation.
 8011 Continue
      Restart=.True.
      Call Qpg_iScalar('LoProp Restart',Found)
      If (Found) Then
         Call Get_iScalar('LoProp Restart',iRestart)
      End If
      If (iRestart .eq. 0) Then
         Write(6,*) 'LoProp was not able to restart.'
         Write(6,*) 'Make sure that LoProp was completed on a previous'
     &            //'run.'
         Call Quit_OnUserError()
      End If
      GoTo 999

*
*>>>>>>>>>>>>>> TDEN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Tell LoProp to collect a transition density from a previous
* Rassi-calculation.
 8012 Continue
      TDensity=.true.
      NoField=.true.
      Key=Get_Ln(LuSpool)
      Call Get_I1(1,nStateI)
      Call Get_I1(2,nStateF)
      Go To 999

*
*>>>>>>>>>>>>>>> XHOLe <<<<<<<<<<<<<<<<<<<<<<<<<<<<
* Compute and distribute exchange-hole dipole moments for
* dispersion coefficients.
 8013 Continue
      XHole=.true.
      NoField=.true.
      Go To 999

*
*>>>>>>>>>>>>>>>> DIFFuse <<<<<<<<<<<<<<<<<<<<<<<<<<
* Section for turning the LoProp moments into diffuse
* functions, i.e. obtain an exponent. No moving of bond
* stuff allowed.
 8014 Continue
      NoField=.true.
      Bond_Threshold=1.0d20
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)
      If(Line(1:4).eq.'NUME') then
        Diffuse(1)=.true.
        Diffuse(2)=.true.
80141   Continue
        Key=Get_Ln(LuSpool)
        Line = Key
        Call UpCase(Line)
        If(Line(1:4).eq.'LIMI') then
          Key = Get_Ln(LuSpool)
          Call Get_F(1,dLimmo,2)
        Elseif(Line(1:4).eq.'THRE') then
          Key=Get_Ln(LuSpool)
          Call Get_F1(1,Thrs1)
          Call Get_F1(2,Thrs2)
          Call Get_I1(3,nThrs)
          Call Get_F1(4,ThrsMul)
        Elseif(Line(1:4).eq.'END ') then
          GoTo 999
        Else
          Write(6,*) 'Undefined option for ''DIFFuse'':',Key
          Call FindErrorLine
          Call Quit_OnUserError()
        Endif
        GoTo 80141
      Elseif(Line(1:4).eq.'REXT') then
        Diffuse(1)=.true.
        Diffuse(3)=.true.
      Else
        Write(6,*) 'Undefined option for ''DIFFuse'':',Key
        Call FindErrorLine
        Call Quit_OnUserError()
      Endif
      Go To 999
*>>>>>>>>>>>>> ALPH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 8015 Continue
* Change the alpha in the penalty function for the
* fluctuating charge contribution to polarisabilities
      Key = Get_Ln(LuSpool)
      Call Get_F1(1,Alpha)
      GoTo 999

 8016 Continue
      LIonize=.true.
      GoTo 999
*
*>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 9000 Continue
*
      Write (6,*)
      If (NoField) Then
         Write (6,*) ' No dynamic properties will be computed.'
         Write (6,*)
      Else
         If (Restart .AND. iRestart .eq. 2) Then
            Write(6,*) ' Previous LoProp calculation was run with the '
     &               //'NOFIeld flag.'
            Write(6,*) ' Thus it is not possible to restart and '
     &               //'calculate dynamic properties.'
            Call Quit_OnUserError()
         End If
         Write (6,*) ' Dynamic properties will be computed.'
         Write (6,*)
         Write (6,'(A,F12.6,A)') '  Applied field +/-', Delta,' au'
         Write (6,*)
      End If
      Write (6,*) ' Expansion centers of the domains are for the'
      If (Standard) Then
         Write (6,*) '  atomic domains: the atom center'
         Write (6,*) '  bond domains  : the center of the bond'
      Else
         Write (6,*) '  atomic domains: the center which set the dipole'
     &               //' moments to zero'
         Write (6,*) '  bond domains  : the center which minimize the '
     &               //'diagonal terms of the quadrupole moment'
         Write (6,*)
         Write (6,*) ' Observe that if the first non-zero term in the '
     &               //'expansion does not dominate the centers are '
     &               //'the original atomic and bond centers!'
      End If
      Write (6,*)
      If(UserDen) then
        Write(6,*)' Read density matrix from user.'
        Write(6,*)
      Endif
      If(TDensity) then
        Write(6,*)' Use transition density matrix from Rassi.'
        Write(6,*)
      Endif
      If(XHole) then
        Write(6,*)' Exchange hole second moment computation and '
     &          //'localization.'
        Write(6,*)
      Endif
      If(Diffuse(1)) then
        Write(6,*)' Computation of exponents to non-zero width '
     &          //'Slater functions.'
        If(Diffuse(2)) then
          Write(6,*)' --- Numerical determination.'
        Elseif(Diffuse(3)) then
          Write(6,*)' --- Analytical determination.'
        Endif
        Write(6,*)
      Endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
