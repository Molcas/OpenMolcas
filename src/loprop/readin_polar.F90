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

subroutine ReadIn_Polar(NoField,Delta,MpProp_Level,Bond_Threshold,iPlot,iPrint,Standard,Opt_Method,UserDen,PrintDen,SubtractDen, &
                        SubScale,Restart,TDensity,nStateI,nStateF,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,Alpha,LIonize)

use Constants, only: One, Two, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(inout) :: NoField
real(kind=wp), intent(out) :: Delta, Bond_Threshold, SubScale, dLimmo(2), Thrs1, Thrs2, ThrsMul, Alpha
integer(kind=iwp), intent(out) :: MpProp_Level, iPlot, iPrint, nStateI, nStateF, nThrs
logical(kind=iwp), intent(out) :: Standard, UserDen, PrintDen, SubtractDen, Restart, TDensity, Diffuse(3), LIonize
character(len=12), intent(out) :: Opt_Method
!---- Define local variables
integer(kind=iwp) :: LuSpool, iRestart, lMax
logical(kind=iwp) :: Found
character(len=180) :: Key, Line
character(len=180), external :: Get_Ln

! copy input from standard input to a local scratch file

LuSpool = 21
call SpoolInp(LuSpool)

! read input

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

Delta = 0.001_wp
call Get_iScalar('Highest Mltpl',lMax)
MpProp_Level = lMax
Bond_Threshold = OneHalf
iPlot = 0
iPrint = 0
Standard = .true.
UserDen = .false.
PrintDen = .false.
SubtractDen = .false.
Restart = .false.
iRestart = 0
SubScale = One
nStateI = 1
nStateF = 2
Opt_Method = ' '
TDensity = .false.
Diffuse(:) = .false.
dLimmo(:) = [0.65_wp,Two]
Thrs1 = 1.0e-5_wp
Thrs2 = 1.0e-4_wp
nThrs = 3
ThrsMul = 1.0e-2_wp
Alpha = 7.1421297_wp
LIonize = .false.

! Comment on Alpha:
! This value of Alpha is for backward-compability.
! For large systems it may have to be reduced, for
! example to 2.0.

!---- Locate "start of input"
rewind(LuSpool)
call RdNLst(LuSpool,'LoProp')

do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)

  select case (Line(1:4))
    !case ('TIT')
      !>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    case ('NOFI')
      !>>>>>>>>>>>>> NOFI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      NoField = .true.

    case ('DELT')
      !>>>>>>>>>>>>> DELT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      Key = Get_Ln(LuSpool)
      call Get_F1(1,Delta)

    case ('EXPA')
      !>>>>>>>>>>>>> EXPA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        write(u6,*) 'Undefined option for "EXPAnsion center":',Key
        call FindErrorLine()
        call Quit_OnUserError()
      end if

    case ('MPPR')
      !>>>>>>>>>>>>> MPPR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Read max multipole level for output in the MpProp file
      Key = Get_Ln(LuSpool)
      call Get_I1(1,MpProp_Level)

    case ('BOND')
      !>>>>>>>>>>>>> BOND <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Read max bond length - all bonds longer than this will be ignored
      Key = Get_Ln(LuSpool)
      call Get_F1(1,Bond_Threshold)

    case ('PLOT')
      !>>>>>>>>>>>>> PLOT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Provide informations needed for plotting t vs. bond coordinate
      ! as well as printing a table of t values.
      iPlot = 1

    case ('PRIN')
      !>>>>>>>>>>>>> PRIN   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Read print level
      Key = Get_Ln(LuSpool)
      call Get_I1(1,iPrint)

    case ('USER')
      !>>>>>>>>>>>>>> USER   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Tell LoProp to read in density matrix(es) supplied by user.
      UserDen = .true.

    case ('PRDE')
      !>>>>>>>>>>>>>> PRDE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Tell LoProp to output density matrix to file
      PrintDen = .true.
      NoField = .true. !Only static properties allowed

    case ('SUBD')
      !>>>>>>>>>>>>>> SUBD   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Tell LoProp to subtract a user supplied density matrix
      ! from the current one (which could also be a USER supplied
      ! density matrix). Also input a scaling factor which the
      ! difference density is multiplied by.
      SubtractDen = .true.
      NoField = .true. !Only static properties allowed
      Key = Get_Ln(LuSpool)
      call Get_F1(1,SubScale)

    case ('REST')
      !>>>>>>>>>>>>> REST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Restart LoProp with the info stored in a previous LoProp
      ! calculation.
      Restart = .true.
      call Qpg_iScalar('LoProp Restart',Found)
      if (Found) then
        call Get_iScalar('LoProp Restart',iRestart)
      end if
      if (iRestart == 0) then
        write(u6,*) 'LoProp was not able to restart.'
        write(u6,*) 'Make sure that LoProp was completed on a previous run.'
        call Quit_OnUserError()
      end if

    case ('TDEN')
      !>>>>>>>>>>>>>> TDEN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Tell LoProp to collect a transition density from a previous
      ! Rassi-calculation.
      TDensity = .true.
      NoField = .true.
      Key = Get_Ln(LuSpool)
      call Get_I1(1,nStateI)
      call Get_I1(2,nStateF)

    case ('DIFF')
      !>>>>>>>>>>>>>>>> DIFFuse <<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Section for turning the LoProp moments into diffuse
      ! functions, i.e. obtain an exponent. No moving of bond
      ! stuff allowed.
      NoField = .true.
      Bond_Threshold = huge(Bond_Threshold)
      Key = Get_Ln(LuSpool)
      Line = Key
      call UpCase(Line)
      if (Line(1:4) == 'NUME') then
        Diffuse(1) = .true.
        Diffuse(2) = .true.
        do
          Key = Get_Ln(LuSpool)
          Line = Key
          call UpCase(Line)
          if (Line(1:4) == 'LIMI') then
            Key = Get_Ln(LuSpool)
            call Get_F(1,dLimmo,2)
          else if (Line(1:4) == 'THRE') then
            Key = Get_Ln(LuSpool)
            call Get_F1(1,Thrs1)
            call Get_F1(2,Thrs2)
            call Get_I1(3,nThrs)
            call Get_F1(4,ThrsMul)
          else if (Line(1:4) == 'END ') then
            exit
          else
            write(u6,*) 'Undefined option for "DIFFuse":',Key
            call FindErrorLine()
            call Quit_OnUserError()
          end if
        end do
      else if (Line(1:4) == 'REXT') then
        Diffuse(1) = .true.
        Diffuse(3) = .true.
      else
        write(u6,*) 'Undefined option for "DIFFuse":',Key
        call FindErrorLine()
        call Quit_OnUserError()
      end if

    case ('ALPH')
      !>>>>>>>>>>>>> ALPH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! Change the alpha in the penalty function for the
      ! fluctuating charge contribution to polarisabilities
      Key = Get_Ln(LuSpool)
      call Get_F1(1,Alpha)

    case ('LION')
      !>>>>>>>>>>>>> LION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      LIonize = .true.

    case ('END ')
      !>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      exit

    case default
      write(u6,*) 'Unidentified key word:',Key
      call FindErrorLine()
      call Quit_OnUserError()

  end select

end do

write(u6,*)
if (NoField) then
  write(u6,*) ' No dynamic properties will be computed.'
  write(u6,*)
else
  if (Restart .and. (iRestart == 2)) then
    write(u6,*) ' Previous LoProp calculation was run with the NOFIeld flag.'
    write(u6,*) ' Thus it is not possible to restart and calculate dynamic properties.'
    call Quit_OnUserError()
  end if
  write(u6,*) ' Dynamic properties will be computed.'
  write(u6,*)
  write(u6,'(A,F12.6,A)') '  Applied field +/-',Delta,' au'
  write(u6,*)
end if
write(u6,*) ' Expansion centers of the domains are for the'
if (Standard) then
  write(u6,*) '  atomic domains: the atom center'
  write(u6,*) '  bond domains  : the center of the bond'
else
  write(u6,*) '  atomic domains: the center which set the dipole moments to zero'
  write(u6,*) '  bond domains  : the center which minimize the diagonal terms of the quadrupole moment'
  write(u6,*)
  write(u6,*) ' Observe that if the first non-zero term in the expansion does not dominate,'
  write(u6,*) ' the centers are the original atomic and bond centers!'
end if
write(u6,*)
if (UserDen) then
  write(u6,*) ' Read density matrix from user.'
  write(u6,*)
end if
if (TDensity) then
  write(u6,*) ' Use transition density matrix from Rassi.'
  write(u6,*)
end if
if (Diffuse(1)) then
  write(u6,*) ' Computation of exponents to non-zero width Slater functions.'
  if (Diffuse(2)) then
    write(u6,*) ' --- Numerical determination.'
  else if (Diffuse(3)) then
    write(u6,*) ' --- Analytical determination.'
  end if
  write(u6,*)
end if

! remove local copy of standard input

call Close_LuSpool(LuSpool)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine ReadIn_Polar
