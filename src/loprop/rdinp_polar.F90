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

subroutine RdInp_Polar(LuSpool,NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen, &
                       SubtractDen,SubScale,Restart,TDensity,nStateI,nStateF,XHole,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul, &
                       Alpha,LIonize)

implicit real*8(a-h,o-z)
!---- Define local variables
character*180 Key, Line
character*180 Get_Ln
character*12 Opt_Method
external Get_Ln
logical NoField, Standard, UserDen, PrintDen, SubtractDen
logical Restart, Found, TDensity, XHole, Diffuse(3)
logical LIonize
dimension dLimmo(2)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

Delta = 0.001d00
call Get_iScalar('Highest Mltpl',lMax)
MpProp_Level = lMax
Bond_Threshold = 1.5d0
iPlot = 0
iPrint = 0
Standard = .true.
UserDen = .false.
PrintDen = .false.
SubtractDen = .false.
Restart = .false.
iRestart = 0
SubScale = 1.0d0
nStateI = 1
nStateF = 2
Opt_Method = ' '
TDensity = .false.
XHole = .false.
Diffuse(1) = .false.
Diffuse(2) = .false.
Diffuse(3) = .false.
dLimmo(1) = 0.65d0
dLimmo(2) = 2.0d0
Thrs1 = 1d-5
Thrs2 = 1d-4
nThrs = 3
ThrsMul = 1d-2
Alpha = 7.1421297d0
LIonize = .false.

! Comment on Alpha:
! This value of Alpha is for backward-compability.
! For large systems it may have to be reduced, for
! example to 2.0.

!---- Locate "start of input"
rewind(LuSpool)
call RdNLst(LuSpool,'LoProp')

999 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)

!if (Line(1:4) == 'TITL') Go To 8000
if (Line(1:4) == 'NOFI') Go To 8001
if (Line(1:4) == 'DELT') Go To 8002
if (Line(1:4) == 'EXPA') Go To 8003
if (Line(1:4) == 'MPPR') Go To 8004
if (Line(1:4) == 'BOND') Go To 8005
if (Line(1:4) == 'PLOT') Go To 8006
if (Line(1:4) == 'PRIN') Go To 8007
if (Line(1:4) == 'USER') Go To 8008
if (Line(1:4) == 'PRDE') Go To 8009
if (Line(1:4) == 'SUBD') Go To 8010
if (Line(1:4) == 'REST') Go To 8011
if (Line(1:4) == 'TDEN') Go To 8012
if (Line(1:4) == 'XHOL') Go To 8013
if (Line(1:4) == 'DIFF') Go To 8014
if (Line(1:4) == 'ALPH') Go To 8015
if (Line(1:4) == 'LION') Go To 8016
if (Line(1:4) == 'END ') Go To 9000
write(6,*) 'Unidentified key word:',Key
call FindErrorLine
call Quit_OnUserError()

!>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!8000 Continue
goto 999

!>>>>>>>>>>>>> NOFI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8001 continue
NoField = .true.
goto 999

!>>>>>>>>>>>>> DELT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8002 continue
Key = Get_Ln(LuSpool)
call Get_F1(1,Delta)
goto 999

!>>>>>>>>>>>>> EXPA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8003 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)

if (Line(1:4) == 'MIDP') then
  Standard = .true.
else if (Line(1:4) == 'OPTI') then
  Standard = .false.
  Opt_Method = 'Optimized'
else if (Line(1:4) == 'MULT') then
  Standard = .false.
  Opt_Method = 'Multipole'
else
  write(6,*) 'Undefined option for ''EXPAnsion center'':',Key
  call FindErrorLine
  call Quit_OnUserError()
end if
goto 999

!>>>>>>>>>>>>> MPPR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Read max multipole level for output in the MpProp file
8004 continue
Key = Get_Ln(LuSpool)
call Get_I1(1,MpProp_Level)
goto 999

!>>>>>>>>>>>>> BOND <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Read max bond length - all bonds longer than this will be ignored
8005 continue
Key = Get_Ln(LuSpool)
call Get_F1(1,Bond_Threshold)
goto 999

!>>>>>>>>>>>>> PLOT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Provide informations needed for plotting t vs. bond coordinate
! as well as printing a table of t values.
8006 continue
iPlot = 1
goto 999

!>>>>>>>>>>>>> PRIN   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Read print level
8007 continue
Key = Get_Ln(LuSpool)
call Get_I1(1,iPrint)
goto 999

!>>>>>>>>>>>>>> USER   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Tell LoProp to read in density matrix(es) supplied by user.
8008 continue
UserDen = .true.
goto 999

!>>>>>>>>>>>>>> PRDE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Tell LoProp to output density matrix to file
8009 continue
PrintDen = .true.
NoField = .true. !Only static properties allowed
goto 999

!>>>>>>>>>>>>>> SUBD   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Tell LoProp to subtract a user supplied density matrix
! from the current one (which could also be a USER supplied
! density matrix). Also input a scaling factor which the
! difference density is multiplied by.
8010 continue
SubtractDen = .true.
NoField = .true. !Only static properties allowed
Key = Get_Ln(LuSpool)
call Get_F1(1,SubScale)
goto 999

!>>>>>>>>>>>>> NOFI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Restart LoProp with the info stored in a previous LoProp
! calculation.
8011 continue
Restart = .true.
call Qpg_iScalar('LoProp Restart',Found)
if (Found) then
  call Get_iScalar('LoProp Restart',iRestart)
end if
if (iRestart == 0) then
  write(6,*) 'LoProp was not able to restart.'
  write(6,*) 'Make sure that LoProp was completed on a previous run.'
  call Quit_OnUserError()
end if
goto 999

!>>>>>>>>>>>>>> TDEN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Tell LoProp to collect a transition density from a previous
! Rassi-calculation.
8012 continue
TDensity = .true.
NoField = .true.
Key = Get_Ln(LuSpool)
call Get_I1(1,nStateI)
call Get_I1(2,nStateF)
Go To 999

!>>>>>>>>>>>>>>> XHOLe <<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Compute and distribute exchange-hole dipole moments for
! dispersion coefficients.
8013 continue
XHole = .true.
NoField = .true.
Go To 999

!>>>>>>>>>>>>>>>> DIFFuse <<<<<<<<<<<<<<<<<<<<<<<<<<
! Section for turning the LoProp moments into diffuse
! functions, i.e. obtain an exponent. No moving of bond
! stuff allowed.
8014 continue
NoField = .true.
Bond_Threshold = 1.0d20
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)
if (Line(1:4) == 'NUME') then
  Diffuse(1) = .true.
  Diffuse(2) = .true.
80141 continue
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  if (Line(1:4) == 'LIMI') then
    Key = Get_Ln(LuSpool)
    call Get_F(1,dLimmo,2)
  elseif (Line(1:4) == 'THRE') then
    Key = Get_Ln(LuSpool)
    call Get_F1(1,Thrs1)
    call Get_F1(2,Thrs2)
    call Get_I1(3,nThrs)
    call Get_F1(4,ThrsMul)
  elseif (Line(1:4) == 'END ') then
    goto 999
  else
    write(6,*) 'Undefined option for ''DIFFuse'':',Key
    call FindErrorLine
    call Quit_OnUserError()
  end if
  goto 80141
elseif (Line(1:4) == 'REXT') then
  Diffuse(1) = .true.
  Diffuse(3) = .true.
else
  write(6,*) 'Undefined option for ''DIFFuse'':',Key
  call FindErrorLine
  call Quit_OnUserError()
end if
Go To 999

!>>>>>>>>>>>>> ALPH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8015 continue
! Change the alpha in the penalty function for the
! fluctuating charge contribution to polarisabilities
Key = Get_Ln(LuSpool)
call Get_F1(1,Alpha)
goto 999

8016 continue
LIonize = .true.
goto 999

!>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
9000 continue

write(6,*)
if (NoField) then
  write(6,*) ' No dynamic properties will be computed.'
  write(6,*)
else
  if (Restart .and. iRestart == 2) then
    write(6,*) ' Previous LoProp calculation was run with the NOFIeld flag.'
    write(6,*) ' Thus it is not possible to restart and calculate dynamic properties.'
    call Quit_OnUserError()
  end if
  write(6,*) ' Dynamic properties will be computed.'
  write(6,*)
  write(6,'(A,F12.6,A)') '  Applied field +/-',Delta,' au'
  write(6,*)
end if
write(6,*) ' Expansion centers of the domains are for the'
if (Standard) then
  write(6,*) '  atomic domains: the atom center'
  write(6,*) '  bond domains  : the center of the bond'
else
  write(6,*) '  atomic domains: the center which set the dipole moments to zero'
  write(6,*) '  bond domains  : the center which minimize the diagonal terms of the quadrupole moment'
  write(6,*)
  write(6,*) ' Observe that if the first non-zero term in the expansion does not dominate,'
  write(6,*) ' the centers are the original atomic and bond centers!'
end if
write(6,*)
if (UserDen) then
  write(6,*) ' Read density matrix from user.'
  write(6,*)
end if
if (TDensity) then
  write(6,*) ' Use transition density matrix from Rassi.'
  write(6,*)
end if
if (XHole) then
  write(6,*) ' Exchange hole second moment computation and localization.'
  write(6,*)
end if
if (Diffuse(1)) then
  write(6,*) ' Computation of exponents to non-zero width Slater functions.'
  if (Diffuse(2)) then
    write(6,*) ' --- Numerical determination.'
  elseif (Diffuse(3)) then
    write(6,*) ' --- Analytical determination.'
  end if
  write(6,*)
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine RdInp_Polar
