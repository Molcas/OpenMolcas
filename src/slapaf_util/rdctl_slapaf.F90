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

subroutine RdCtl_Slapaf(LuSpool,Dummy_Call)

use kriging_mod
use ThermoChem
use Symmetry_Info, only: Symmetry_Info_Get
use Slapaf_Info, only: Cx, Gx, Weights, MF, Atom, nSup, RefGeo, GradRef, nStab, Lbl, mRowH, Coor
use Slapaf_Parameters, only: iRow, iRow_c, ddV_Schlegel, HWRS, iOptH, HrmFrq_Show, IRC, Curvilinear, Redundant, FindTS, nBVec, &
                             User_Def, MaxItr, iOptC, rHidden, CnstWght, lOld, Beta, Beta_Disp, Line_Search, TSConstraints, &
                             GNrm_Threshold, Mode, ThrEne, ThrGrd, nLambda, ThrCons, ThrMEP, Baker, eMEPTest, rMEP, MEP, nMEP, &
                             MEPNum, MEPCons, dMEPStep, MEP_Type, MEP_Algo, Max_Center, Delta, RtRnc, rFuzz, lNmHss, Cubic, &
                             Request_Alaska, CallLast, lCtoF, Track, isFalcon, MxItr, nWndw, Iter, WeightedConstraints, NADC, &
                             Fallback
use UnixInfo, only: SuperName

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
integer iDum(1)
logical Found, Dummy_Call
character(len=180) Get_Ln
character*16 FilNam
character*3 MEPLab
character(len=180), parameter :: BLine = ''
character(len=180) :: Key = '', Char = ''
real*8, allocatable :: DIR(:,:), Tmp(:), TmpRx(:)
#include "cgetl.fh"
external Get_Ln
logical External_UDC, Explicit_IRC, Expert, ThrInp, FirstNum, Manual_Beta
#include "angstr.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
Expert = .false.
Lu = 6
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate some parameters

call Symmetry_Info_Get()
call Init_Slapaf()
nsAtom = size(Coor,2)
iPrint = nPrint(iRout)
iSetAll = 2**30-1

call f_Inquire('UDC.Gateway',External_UDC)

iMEP = 0
Explicit_IRC = .false.
WeightedConstraints = .false.
ThrInp = .false.
call Qpg_iScalar('nMEP',Found)
if (Found) call Get_iScalar('nMEP',iMEP)
if (iMEP == 0) then
  iOff_Iter = 0
  call Put_iScalar('iOff_Iter',iOff_Iter)
end if
Manual_Beta = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! When called from outside Slapaf or as a dummy call, process no
! input but proceed with default values only.

if ((SuperName /= 'slapaf') .or. Dummy_Call) then
  Char = 'END '
  Go To 666
end if
!                                                                      *
!***********************************************************************
!*************************   Input section   ***************************
!***********************************************************************
!                                                                      *
LuRd = LuSpool
call RdNlst(LuRd,'SLAPAF')
999 Char = Get_Ln(LuRd)
666 continue
call UpCase(Char)
!write(Lu,'(A)') Char
!write(Lu,*) iOptC
if (Char == BLine) Go To 999
if (char(1:1) == '*') Go To 999
!if (Char(1:4) == 'AIL ') Go To 102
!if (Char(1:4) == 'AIP ') Go To 104
!if (Char(1:4) == 'AISP') Go To 105
!if (Char(1:4) == 'AIME') Go To 107
!if (Char(1:4) == 'AIBL') Go To 108
!if (Char(1:4) == 'AIMB') Go To 109
!if (Char(1:4) == 'L-VA') Go To 112
if (char(1:4) == 'NDEL') Go To 113
if (char(1:4) == 'BAKE') Go To 926
if (char(1:4) == 'C1-D') Go To 936
if (char(1:4) == 'C2-D') Go To 937
if (char(1:4) == 'CART') Go To 918
if (char(1:4) == 'CNWE') Go To 990
if (char(1:4) == 'CONS') Go To 9478
if (char(1:4) == 'CTOF') Go To 904
if (char(1:4) == 'CUBI') Go To 947
if (char(1:4) == 'DDVS') Go To 9271
if (char(1:4) == 'DELT') Go To 946
if (char(1:4) == 'DISO') Go To 9452
if (char(1:4) == 'DXDX') Go To 939
if (char(1:4) == 'DXG ') Go To 940
if (char(1:4) == 'END ') Go To 998
if (char(1:4) == 'EXPE') Go To 993
if (char(1:4) == 'EXTR') Go To 971
if (char(1:4) == 'FALC') Go To 800
if (char(1:4) == 'FIND') Go To 963
if (char(1:4) == 'FUZZ') Go To 123
if (char(1:4) == 'GDX ') Go To 940
if (char(1:4) == 'GG  ') Go To 941
if (char(1:4) == 'GNRM') Go To 968
if (char(1:4) == 'GRAD') Go To 979
if (char(1:4) == 'HRMS') Go To 995
if (char(1:4) == 'HUPD') Go To 914
if (char(1:4) == 'HWRS') Go To 929
if (char(1:4) == 'INTE') Go To 902
if (char(1:4) == 'IRC ') Go To 997
if (char(1:4) == 'ITER') Go To 925
if (char(1:4) == 'KRIG') Go To 100
if (char(1:4) == 'LAST') Go To 9280
if (char(1:4) == 'LINE') Go To 9281
if (char(1:4) == 'MAXS') Go To 915
if (char(1:4) == 'MAXD') Go To 916
if ((char(1:4) == 'MEP-') .or. (char(1:4) == 'MEP ')) Go To 964
if ((char(1:4) == 'MEPA') .or. (char(1:4) == 'IRCA')) Go To 322
if ((char(1:4) == 'MEPC') .or. (char(1:4) == 'IRCC')) Go To 323
if ((char(1:4) == 'MEPS') .or. (char(1:4) == 'IRCS')) Go To 9971
if ((char(1:4) == 'MEPT') .or. (char(1:4) == 'IRCT')) Go To 321
if (char(1:4) == 'MODE') Go To 942
if (char(1:4) == 'MXMI') Go To 106
if ((char(1:4) == 'NMEP') .or. (char(1:4) == 'NIRC')) Go To 965
if (char(1:4) == 'NEWT') Go To 935
if (char(1:4) == 'NOEM') Go To 991
if (char(1:4) == 'NOFA') Go To 114
if (char(1:4) == 'NOHW') Go To 960
if (char(1:4) == 'NOLA') Go To 930
if (char(1:4) == 'NOLI') Go To 928
if (char(1:4) == 'NOWB') Go To 984
if (char(1:4) == 'NOWC') Go To 985
if (char(1:4) == 'NOWH') Go To 986
if (char(1:4) == 'NUME') Go To 945
if (char(1:4) == 'OLDF') Go To 903
if (char(1:4) == 'PRFC') Go To 9201
if (char(1:4) == 'PRIN') Go To 920
if (char(1:4) == 'RATI') Go To 938
if (char(1:4) == 'REDU') Go To 994
if (char(1:4) == 'REAC') Go To 996
if (char(1:4) == 'REFE') Go To 966
if (char(1:4) == 'RHID') Go To 988
if (char(1:4) == 'RMEP') Go To 980
if (char(1:4) == 'RS-P') Go To 967
if (char(1:4) == 'RTRN') Go To 962
if (char(1:4) == 'SUPS') Go To 911
if (char(1:4) == 'TFOF') Go To 110
if (char(1:4) == 'THER') Go To 9451
if (char(1:4) == 'THRS') Go To 908
if (char(1:4) == 'TOLE') Go To 909
if (char(1:4) == 'TRAC') Go To 910
if (char(1:4) == 'TS  ') Go To 951
if (char(1:4) == 'TSCO') Go To 320
if (char(1:4) == 'VDWB') Go To 981
if (char(1:4) == 'VDWC') Go To 982
if (char(1:4) == 'VDWH') Go To 983
if (char(1:4) == 'WIND') Go To 934
call WarningMessage(2,'Error in RdCtl_Slapaf')
if (char(1:1) == ' ') then
  write(Lu,*) ' RdCtl_Slapaf: Command line starts with a blank.'
else
  write(Lu,*)
  write(Lu,*) ' *********** ERROR ***********'
  write(Lu,*) ' The program has been supplied'
  write(Lu,*) ' with an unknown command.     '
  write(Lu,*) ' *****************************'
end if
write(Lu,'(A)') Char
call Quit_OnUserError()
!                                                                      *
!***** INTE ************************************************************
!                                                                      *
! Read the internal coordinate specification.

902 continue
New_Line = 1
Lu_UDIC = 91
FilNam = 'UDIC'
call molcas_open(Lu_UDIC,FilNam)
rewind(Lu_UDIC)

! mInt is the number of internal coordinates you will define.
! mInt = nDimBC - mTROld
! Subroutine DefInt defines the B matrix.
! The matrix B relates a shift in an internal coordinate to
! shifts in cartesian coordinates,
!
!           |dq> = B |dx>
!                      =
! and has the dimension (3*nsAtom x mInt).
992 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
if (Key(1:4) == 'END ') then
  close(Lu_UDIC)
  Go To 999
end if

! Here is a fix because auto will break up the lines if there is an
! equal sign in the input.

! Lines with VARY or FIX doesn't have equal signs

if (Key(1:4) == 'VARY') nBVec = iRow
if ((Key(1:4) == 'VARY') .or. (Key(1:3) == 'FIX') .or. (Key(1:4) == 'ROWH')) then
  New_Line = 0
end if

111 continue
if (New_Line == 1) then
  if (index(Key,'=') == 0) call FixEqualSign2(Key,LuRd,Lu_UDIC,iRow,New_Line)
  if (New_Line == 2) then
    close(Lu_UDIC)
    Go To 999
  end if
  Go To 111
end if

iRow = iRow+1

write(Lu_UDIC,'(A)') Key

! If this line does not have a continuation the next line should
! have a equal sign!
if (index(Key,'&') == 0) New_Line = 1
Go To 992
!                                                                      *
!***** CTOF ************************************************************
!                                                                      *
! Read the internal (C)oordinate specification (TO) be (F)ollowed.

904 continue
if (iRow > 0) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '*********** ERROR ***************'
  write(Lu,*) 'CtoF and User-defined Coordinates'
  write(Lu,*) 'are mutually exclusive.          '
  write(Lu,*) '*********************************'
  call Quit_OnUserError()
end if
iNull = 0
New_Line = 1
lCtoF = .true.
Lu_UDIC = 91
FilNam = 'UDIC'
call molcas_open(Lu_UDIC,FilNam)
rewind(Lu_UDIC)
Key = Get_Ln(LuRd)
call UpCase(Key)
call FixEqualSign2(Key,LuRd,Lu_UDIC,iNull,New_Line)
write(Lu_UDIC,'(A)') Key
close(Lu_UDIC)
Go To 999
!                                                                      *
!***** CONS ************************************************************
!                                                                      *
! Copy constraints definition into the UDC file, to be read
! (after fixing and merging, if necessary) by DefInt2/Cllct2.

9478 continue
if (.not. Expert) then
  write(Lu,*)
  write(Lu,*) ' ************ ERROR ***************'
  write(Lu,*) ' Obsolete input standard!'
  write(Lu,*) ' The CONSTRAINT section should'
  write(Lu,*) ' be define in the &Gateway input.'
  write(Lu,*)
  write(Lu,*) ' To override add the EXPERT keyword'
  write(Lu,*) ' to the top of the SLAPAF input.'
  write(Lu,*) ' **********************************'
  call Quit_OnUserError()
end if
if (External_UDC) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '****************** ERROR *********************'
  write(Lu,*) 'Multiple definitions of constraints.'
  write(Lu,*) 'Check Gateway and Slapaf inputs for conflicts!'
  write(Lu,*) '**********************************************'
  call Quit_OnUserError()
end if
Lu_UDC = 20
FilNam = 'UDC'
Lu_UDC = IsFreeUnit(Lu_UDC)
call Molcas_Open(Lu_UDC,FilNam)
318 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
Key = adjustl(Key)
write(Lu_UDC,'(A)') trim(Key)
if (Key(1:4) /= 'END') Go To 318
close(Lu_UDC)
Go To 999
!                                                                      *
!***** VDWB VdW correction both coordinate and Hessian *****************
!                                                                      *
981 iOptC = ior(1024,iOptC)
iOptC = ior(2048,iOptC)
Go To 999
!                                                                      *
!***** NO VDWB VdW correction both coordinate and Hessian **************
!                                                                      *
984 Mask = iSetAll-2**10-2**11
iOptC = iand(Mask,iOptC)
Go To 999
!                                                                      *
!***** VDWB VdW correction for coordinate only *************************
!                                                                      *
982 iOptC = ior(2048,iOptC)
Go To 999
!                                                                      *
!***** NO VDWB VdW correction for coordinate only **********************
!                                                                      *
985 Mask = iSetAll-2**11
iOptC = iand(Mask,iOptC)
Go To 999
!                                                                      *
!***** VDWB VdW correction for Hessian only ****************************
!                                                                      *
983 iOptC = ior(1024,iOptC)
Go To 999
!                                                                      *
!***** NO VDWB VdW correction for Hessian only *************************
!                                                                      *
986 Mask = iSetAll-2**10
iOptC = iand(Mask,iOptC)
Go To 999
!                                                                      *
!***** OLDF ************************************************************
!                                                                      *
903 lOld = .true.
Go To 999
!                                                                      *
!***** CART ************************************************************
!                                                                      *
918 CurviLinear = .false.
Go To 999
!                                                                      *
!***** THRS ************************************************************
!                                                                      *
! read the gradient threshold

908 Char = Get_Ln(LuRd)
call Get_F1(1,ThrEne)
call Get_F1(2,ThrGrd)
ThrInp = .true.
Go To 999
!                                                                      *
!***** TOLE ************************************************************
!                                                                      *
! read the constraints threshold

909 Char = Get_Ln(LuRd)
call Get_F1(1,ThrCons)
ThrCons = abs(ThrCons)
Go To 999
!                                                                      *
!***** SUPS ************************************************************
!                                                                      *
! Introduce supersymmetry
! Input format
! nsg                (number of super groups)
! Repeat nsg times
! nmem, (ind.., i = 1, nmem)

911 Char = Get_Ln(LuRd)
call Get_I1(1,nSupSy)
call mma_allocate(nSup,NSUPSY,Label='nSup')
call mma_allocate(Atom,nsAtom,Label='Atom')
jStrt = 1
do i=1,nSupSy
  read(LuRd,*,Err=9630) nSup(i),(Atom(j),j=jStrt,jStrt+nSup(i)-1)
  jStrt = jStrt+nSup(i)
end do
Go To 999
9630 call WarningMessage(2,'Error in RdCtl_Slapaf')
write(Lu,*)
write(Lu,*) '************ ERROR ****************'
write(Lu,*) 'Error while reading supersymmetry.'
write(Lu,*) '***********************************'
call Quit_OnUserError()
!                                                                      *
!***** HUPD ************************************************************
!                                                                      *
914 Char = Get_Ln(LuRd)
read(Char,*) Char
call UpCase(Char)
if (trim(Char) == 'BFGS') then
  iOptH = 4
!else if (trim(Char) == 'MEYER') then
!  iOptH = ior(1,iAnd(iOptH,32))
!else if (trim(Char) == 'BP') then
!  iOptH = ior(2,iAnd(iOptH,32))
else if (trim(Char) == 'NONE') then
  iOptH = ior(8,iand(iOptH,32))
else if (trim(Char) == 'MSP') then
  iOptH = ior(16,iand(iOptH,32))
else if (trim(Char) == 'EU') then
  iOptH = ior(64,iand(iOptH,32))
else if (trim(Char) == 'TS-BFGS') then
  iOptH = ior(128,iand(iOptH,32))
else
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '************ ERROR ****************'
  write(Lu,*) 'Unsupported Hessian update method: ',trim(Char)
  write(Lu,*) '***********************************'
  call Quit_OnUserError()
end if
Go To 999
!                                                                      *
!***** MAXS ************************************************************
!                                                                      *
915 Char = Get_Ln(LuRd)
if (Char == BLine) Go To 915
if (char(1:1) == '*') Go To 915
call Get_F1(1,Beta)
Go To 999
!                                                                      *
!***** MAXD ************************************************************
!                                                                      *
916 Char = Get_Ln(LuRd)
if (Char == BLine) Go To 916
if (char(1:1) == '*') Go To 916
call Get_F1(1,Beta_Disp)
Manual_Beta = .true.
Go To 999
!                                                                      *
!***** PRIN ************************************************************
!                                                                      *
920 Char = Get_Ln(LuRd)
call UpCase(Char)
if (Char == BLine) Go To 920
if (char(1:1) == '*') Go To 920
call Get_I1(1,mPrint)
do i=1,mPrint
922 Char = Get_Ln(LuRd)
  call UpCase(Char)
  if (Char == BLine) Go To 922
  if (char(1:1) == '*') Go To 922
  call Get_I1(1,iRout)
  call Get_I1(2,kPrint)
  nPrint(iRout) = kPrint
end do
Go To 999
!                                                                      *
!***** PRFC ************************************************************
!                                                                      *
! set nPrint to print internal coordinates and hessian

9201 nPrint(21) = 6  ! Eigen-/Values/Vectors of the Hessian (diagmtrx)
nPrint(116) = 6 ! Internal Forces (rlxctl)
if (.not. Request_Alaska) nPrint(30) = 6 ! Coords & Forces (defint)
nPrint(122) = 6 ! Auto-Defined Internal coordinates (printq_sl)
Go To 999
!                                                                      *
!***** ITER ************************************************************
!                                                                      *
! read max iterations

925 Char = Get_Ln(LuRd)
call Get_I1(1,iTmp)
MxItr = min(iTmp,MxItr)
Go To 999
!                                                                      *
!***** KRIG ************************************************************
!                                                                      *
! Activate Kriging

100 Kriging = .true.
Line_Search = .false.
Go To 999
!                                                                      *
!***** NOFA ************************************************************
!                                                                      *
! Deactivate fallback to conventional

114 Fallback = .false.
Go To 999
!                                                                      *
!***** AIMD ************************************************************
!                                                                      *
! Analitical or numerical Matern derivatives
!
!101 Char = Get_Ln(LuRd)
!if ((Char == 'False') .or. (Char == 'false')) then
!  anMd = .false.
!else
!  anMd = .true.
!end if
!Go To 999
!                                                                      *
!***** AIL  ************************************************************
!                                                                      *
! Width limits of the Matern function
!
!102 Char = Get_Ln(LuRd)
!call Get_F(1,lb,3)
!Go To 999
!                                                                      *
!***** AIP  ************************************************************
!                                                                      *
! Parameter of differentiability for Matern function
!
!104 Char = Get_Ln(LuRd)
!call Get_F1(1,pAI)
!if ((pAI > 3) .or. (pAI < 1)) anMd = .false.
!Go To 999
!                                                                      *
!***** AISP ************************************************************
!                                                                      *
! Defining the number of source points for the AI method
!
!105 Char = Get_Ln(LuRd)
!call Get_I1(1,nspAI)
!Go To 999
!                                                                      *
!***** MXMI ************************************************************
!                                                                      *
! Maximum number of Iterations for the Kriging method

106 Char = Get_Ln(LuRd)
call Get_I1(1,Max_Microiterations)
Go To 999
!                                                                      *
!***** AIME ************************************************************
!                                                                      *
! Minimum energy differences of the last two iterations
! (loop exit condition)
!
!107 Char = Get_Ln(LuRd)
!call Get_F1(1,Thr_microiterations)
!Go To 999
!                                                                      *
!***** AIBL ************************************************************
!                                                                      *
! Base line modification value to not ordinary
! (Trend Function on GEK)
!
!108 Char = Get_Ln(LuRd)
!call Get_F1(1,blvAI)
!blAI = .true.
!Go To 999
!                                                                      *
!***** AIMB ************************************************************
!                                                                      *
! Base line modification value to maximum value of the Energy
! This option supersedes any value assigned to blAI
!
!109 Char = Get_Ln(LuRd)
!mblAI = .true.
!Go To 999
!                                                                      *
!***** TFOF ************************************************************
!                                                                      *
! adding energy to the last energy value of the base line
! This option supersedes any value assigned to blAI and mblAI

110 Char = Get_Ln(LuRd)
call Get_F1(1,blavAI)
Go To 999
!                                                                      *
!***** L-VA ************************************************************
!                                                                      *
! Change the l value of the GEK.
!
!112 Char = Get_Ln(LuRd)
!Set_l = .true.
!call Get_F1(1,Value_l)
!call Qpg_dScalar('Value_l',Found)
!if (.not. Found) call Put_dScalar('Value_l',Value_l)
!Go To 999
!                                                                      *
!***** NDELta **********************************************************
!                                                                      *
113 Char = Get_Ln(LuRd)
call Get_I1(1,nD_In)
Go To 999
!                                                                      *
!***** BAKE ************************************************************
!                                                                      *
926 Baker = .true.
Go To 999
!                                                                      *
!***** DDVS ************************************************************
!                                                                      *
9271 DDV_Schlegel = .true.
Go To 999
!                                                                      *
!***** NOLA ************************************************************
!                                                                      *
930 CallLast = .false.
Go To 999
!                                                                      *
!***** NOLI ************************************************************
!                                                                      *
928 Line_Search = .false.
Go To 999
!                                                                      *
!***** LAST ************************************************************
!                                                                      *
9280 Char = Get_Ln(LuRd)
Char = adjustl(Char)
if (Char == BLine) Go To 9280
if (char(1:1) == '*') Go To 9280
call UpCase(Char)
call Put_cArray('LastEnergyMethod',Char,8)
Go To 999
!                                                                      *
!***** LINE ************************************************************
!                                                                      *
9281 Line_Search = .true.
Go To 999
!                                                                      *
!***** HWRS ************************************************************
!                                                                      *
929 HWRS = .true.
Go To 999
!                                                                      *
!***** WIND ************************************************************
!                                                                      *
934 Char = Get_Ln(LuRd)
call UpCase(Char)
if (Char == BLine) Go To 934
if (char(1:1) == '*') Go To 934
call Get_I1(1,nWndw)
Go To 999
!                                                                      *
!***** NEWT ************************************************************
!                                                                      *
935 Mask = iSetAll
Mask = Mask-2**0-2**1-2**2-2**3
iOptC = ior(2**0,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** C1-D ************************************************************
!                                                                      *
936 Mask = iSetAll
Mask = Mask-2**0-2**1-2**2-2**3
iOptC = ior(2**1,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** C2-D ************************************************************
!                                                                      *
937 Mask = iSetAll
Mask = Mask-2**0-2**1-2**2-2**3
iOptC = ior(2**2,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** RATI ************************************************************
!                                                                      *
938 Mask = iSetAll
Mask = Mask-2**0-2**1-2**2-2**3
iOptC = ior(2**3,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** DXDX ************************************************************
!                                                                      *
939 Mask = iSetAll
Mask = Mask-2**4-2**5-2**6
iOptC = ior(2**4,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** DXG  ************************************************************
!                                                                      *
940 Mask = iSetAll
Mask = Mask-2**4-2**5-2**6
iOptC = ior(2**5,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** GG   ************************************************************
!                                                                      *
941 Mask = iSetAll
Mask = Mask-2**4-2**5-2**6
iOptC = ior(2**6,iand(iOptC,Mask))
Go To 999
!                                                                      *
!***** MODE ************************************************************
!                                                                      *
! Mode following algorithm
942 continue
!read(5,'(A)',End=9610) Char
Char = Get_Ln(LuRd)
if (Char == BLine) Go To 942
if (char(1:1) == '*') Go To 942
call Get_I1(1,mode)
Go To 999
!                                                                      *
!***** NUME ************************************************************
!                                                                      *
945 lNmHss = .true.
Go To 999
!                                                                      *
!***** THER ************************************************************
!                                                                      *
9451 lNmHss = .true.
lTherm = .true.
Char = Get_Ln(LuRd)
call Get_I1(1,nsRot)
Char = Get_Ln(LuRd)
call Get_F1(1,UserP)
9454 Char = Get_Ln(LuRd)
call UpCase(Char)
if (char(1:4) == 'END ') then
  if (nUserPT == 0) then
    nUserPT = 1
    UserT(1) = 298.15d0
  end if
  Go To 999
end if
nUserPT = nUserPT+1
call Get_F1(1,UserT(nUserPT))
Go To 9454
!                                                                      *
!***** DISO ************************************************************
!                                                                      *
9452 lDoubleIso = .true.
Go To 999
!                                                                      *
!***** CUBI ************************************************************
!                                                                      *
947 Cubic = .true.
Go To 999
!                                                                      *
!***** DELT ************************************************************
!                                                                      *
946 Char = Get_Ln(LuRd)
call Get_F1(1,Delta)
Go To 999
!                                                                      *
!***** TS   ************************************************************
!                                                                      *
951 Mask = iSetAll-2**7
iOptC = iand(Mask,iOptC)
Go To 999
!                                                                      *
!***** EXTR ************************************************************
!                                                                      *
! Put the program name and the time stamp onto the extract file

971 write(Lu,*) 'RdCtl_Slapaf: EXTRACT option is redundant and is ignored!'
Go To 999
!                                                                      *
!***** NOHW ************************************************************
!                                                                      *
960 HWRS = .false.
Go To 999
!                                                                      *
!***** RTRN ************************************************************
!                                                                      *
962 Char = Get_Ln(LuRd)
call UpCase(Char)
call Get_I1(1,Max_Center)
call Get_F1(2,rtrnc)
if (index(Char,'ANGSTROM') /= 0) Rtrnc = Rtrnc/angstr
Go To 999
!                                                                      *
!***** FIND ************************************************************
!                                                                      *
963 FindTS = .true.
Go To 999
!                                                                      *
!***** TSCO ************************************************************
!                                                                      *
320 LuTS = 20
FilNam = 'TSC'
LuTS = IsFreeUnit(LuTS)
call Molcas_Open(LuTS,FilNam)
319 Key = Get_Ln(LuRd)
call UpCase(Key)
Key = adjustl(Key)
write(LuTS,'(A)') trim(Key)
if (Key(1:4) /= 'END') Go To 319
close(LuTS)
TSConstraints = .true.
Go To 999
!                                                                      *
!***** FUZZ ************************************************************
!                                                                      *
123 Char = Get_Ln(LuRd)
call UpCase(Char)
call Get_F1(1,rFuzz)
if (index(Char,'ANGSTROM') /= 0) rFuzz = rFuzz/angstr
rFuzz = max(rFuzz,1.0D-3)
Go To 999
!                                                                      *
!***** MEP-/MEP  *******************************************************
!                                                                      *
964 MEP = .true.
rMEP = .false.
Go To 999
!                                                                      *
!***** NMEP/NIRC *******************************************************
!                                                                      *
965 Char = Get_Ln(LuRd)
call Get_I1(1,nMEP)
nMEP = min(max(nMEP,1),MaxItr)
Go To 999
!                                                                      *
!***** MEPT/IRCT *******************************************************
!                                                                      *
321 Char = Get_Ln(LuRd)
call UpCase(Char)
if (char(1:6) == 'SPHERE') then
  MEP_Type = 'SPHERE'
else if (char(1:5) == 'PLANE') then
  MEP_Type = 'TRANSVERSE'
else
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '********** ERROR **********'
  write(Lu,*) ' Unrecognized MEP/IRC type.'
  write(Lu,*) '***************************'
  call Quit_OnUserError()
end if
Go To 999
!                                                                      *
!***** MEPA/IRCA *******************************************************
!                                                                      *
322 Char = Get_Ln(LuRd)
call UpCase(Char)
if (char(1:2) == 'GS') then
  MEP_Algo = 'GS'
else if (char(1:2) == 'MB') then
  MEP_Algo = 'MB'
else
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '************* ERROR ************'
  write(Lu,*) ' Unrecognized MEP/IRC algorithm.'
  write(Lu,*) '********************************'
  call Quit_OnUserError()
end if
Go To 999
!                                                                      *
!***** MEPC/IRCC *******************************************************
!                                                                      *
323 Char = Get_Ln(LuRd)
call Get_F1(1,ThrMEP)
ThrMEP = max(Zero,ThrMEP)
Go To 999
!                                                                      *
!***** REFE ************************************************************
!                                                                      *
966 call mma_allocate(RefGeo,3,nsAtom,Label='RefGeo')
call Read_v(LuRd,RefGeo,1,3*nsAtom,1,iErr)
if (iErr /= 0) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '************ ERROR ***************'
  write(Lu,*) 'Error reading reference structure.'
  write(Lu,*) '**********************************'
  call Quit_OnUserError()
end if
Go To 999
!                                                                      *
!***** RS-P ************************************************************
!                                                                      *
967 Mask = iSetAll
Mask = Mask-2**9
iOptC = iand(iOptC,Mask)
Go To 999
!                                                                      *
!***** GNRM ************************************************************
!                                                                      *
968 Char = Get_Ln(LuRd)
call Get_F1(1,GNrm_Threshold)
Go To 999
!                                                                      *
!***** GRAD ************************************************************
!                                                                      *
979 call mma_allocate(GradRef,3,nsAtom,Label='GradRef')

call Read_v(LuRd,GradRef,1,3*nsAtom,1,iErr)
if (iErr /= 0) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '*************** ERROR ******************'
  write(Lu,*) 'Error reading reference gradient vector.'
  write(Lu,*) '****************************************'
  call Quit_OnUserError()
end if

! If there is a transverse vector stored, we are not using this one

call qpg_dArray('Transverse',Found,nRP)
if (Found) call mma_deallocate(GradRef)
Go To 999
!                                                                      *
!***** rMEP ************************************************************
!                                                                      *
980 rMEP = .true.
MEP = .false.
Go To 999
!                                                                      *
!***** rHidden *********************************************************
!                                                                      *
988 Key = Get_Ln(LuRd)
call UpCase(Key)
call Get_F1(1,rHidden)
if (rHidden < Zero) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '************ ERROR *****************'
  write(Lu,*) 'Error reading rHidden. Should be >0.'
  write(Lu,*) '************************************'
  call Quit_OnUserError()
end if
if (index(Key,'ANGSTROM') /= 0) rHidden = rHidden/angstr
Go To 999
!                                                                      *
!***** IRC *************************************************************
!                                                                      *
997 call Qpg_iScalar('IRC',Found)
if (Found) then
  call Get_iScalar('IRC',IRC)
else
  IRC = 1
  call Put_iScalar('IRC',IRC)
end if
MEP = .true.
rMEP = .false.
Go To 999
!                                                                      *
!***** MEPStep/IRCStep *************************************************
!                                                                      *
9971 Char = Get_Ln(LuRd)
call UpCase(Char)
call Get_F1(1,dMEPStep)

! Note that according to the Gonzalez-Schlegel method, only half
! this step is used in the constraint

if (index(Char,'ANGSTROM') /= 0) dMEPStep = dMEPStep/angstr
Go To 999
!                                                                      *
!***** REAC ************************************************************
!                                                                      *
996 Explicit_IRC = .true.
call mma_allocate(TmpRx,3*nsAtom,Label='TmpRx')
call Read_v(LuRd,TmpRx,1,3*nsAtom,1,iErr)
if (IErr /= 0) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '********** ERROR ***********************'
  write(Lu,*) ' Error while reading the Reaction vector'
  write(Lu,*) '****************************************'
  call Quit_OnUserError()
end if
Go To 999
!                                                                      *
!***** HRMS ************************************************************
!                                                                      *
995 HrmFrq_Show = .true.
Go To 999
!                                                                      *
!***** CNWE ************************************************************
!                                                                      *
990 Char = Get_Ln(LuRd)
call Get_F1(1,CnstWght)
Go To 999
!                                                                      *
!***** NOEM ************************************************************
!                                                                      *
991 eMEPTest = .false.
Go To 999
!                                                                      *
!***** EXPE ************************************************************
!                                                                      *
993 Expert = .true.
Go To 999
!                                                                      *
!***** REDU ************************************************************
!                                                                      *
994 Redundant = .true.
Go To 999
!                                                                      *
!***** FALC ************************************************************
!                                                                      *
800 isFalcon = .true.
Go To 999
!                                                                      *
!***** TRAC ************************************************************
!                                                                      *
910 Track = .true.
Go To 999
!                                                                      *
!***********************************************************************
!***********************   End of input section   **********************
!***********************************************************************
!                                                                      *
998 continue
!                                                                      *
!***********************************************************************
!                                                                      *
! Now start fixing constraints.
! First get the external constraints.

if (External_UDC) then
  call Merge_Constraints('UDC.Gateway','','UDC',nLambda,iRow_c)
else
  call Merge_Constraints('','UDC','UDC',nLambda,iRow_c)
end if

! Initial preprocessing

if (iRow_c > 1) then
  Lu_UDC = IsFreeUnit(20)
  call Molcas_Open(Lu_UDC,'UDC')
  call Preprocess_UDC(Lu_UDC,iPrint)
  close(Lu_UDC)
else
  NADC = .false.
end if

! Add NAC if needed

if (NADC) then
  Lu_UDCTMP = IsFreeUnit(20)
  call Molcas_Open(Lu_UDCTMP,'UDCTMP')
  write(Lu_UDCTMP,*) 'NADC = NAC'
  write(Lu_UDCTMP,*) 'VALUE'
  write(Lu_UDCTMP,*) 'NADC = 0.0'
  write(Lu_UDCTMP,*) 'END'
  close(Lu_UDCTMP)
  call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
end if

! Add MEP/IRC if needed

if (MEP .or. rMEP .or. (abs(IRC) == 1)) then
  if (abs(IRC) == 1) then
    MEPLab = 'IRC'
  else
    MEPLab = 'MEP'
  end if
  if (MEPCons .and. (.not. Expert)) then
    call WarningMessage(2,'Error in RdCtl_Slapaf')
    write(Lu,*)
    write(Lu,*) '***************** ERROR ********************'
    write(Lu,*) ' There is a '//trim(Mep_Type)//' constraint that may'
    write(Lu,*) ' conflict with '//MEPLab//' calculations.'
    write(Lu,*) ' You should not explictly specify this constraint,'
    write(Lu,*) ' but just rely on '//MEPLab//'Step/'//MEPLab//'Type keywords.'
    write(Lu,*) ' If you really know what you are doing, you'
    write(Lu,*) ' can use the EXPERT keyword.'
    write(Lu,*) '********************************************'
    call Quit_OnUserError()
  end if
  WeightedConstraints = .true.
  Valu = dMEPStep
  if (MEP_Type == 'SPHERE') Valu = abs(Valu)
  if (MEP .and. (MEP_Algo == 'GS')) Valu = Half*Valu
  if (rMEP) Valu = max(dble(iMEP+1),One)*Valu
  Lu_UDCTMP = IsFreeUnit(20)
  call Molcas_Open(Lu_UDCTMP,'UDCTMP')
  write(Lu_UDCTMP,*) MEPLab//' = '//MEP_Type
  write(Lu_UDCTMP,*) 'VALUE'
  write(Lu_UDCTMP,*) MEPLab//' = ',Valu
  write(Lu_UDCTMP,*) 'END'
  close(Lu_UDCTMP)
  call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
  Beta = min(Beta,abs(Valu))
  MEPnum = nLambda
end if

! Final fixes

call Fix_UDC(iRow_c,nLambda,nsAtom,nStab,.true.)
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate some variables which can only be set after the input has
! been read.

if ((.not. ThrInp) .and. (.not. Baker)) ThrEne = Zero

call Init2()

! Gradients are not needed at the first iteration of a numerical
! Hessian procedure (and only that, i.e. MxItr=0)
FirstNum = (allocated(mRowH) .or. lNmHss .or. Cubic) .and. (Iter == 1) .and. (MxItr == 0)

if ((SuperName == 'slapaf') .and. (.not. FirstNum)) then

  if (Track) then
    call Process_Track()
  else
    call Put_iArray('Root Mapping',iDum,0)
  end if

  call Process_Gradients()

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put in the "Reaction vector" in Cartesians.
! Priority order:
! 1) Explicit by user input (REAC keyword)
! 2) Found on RunOld
! 3) Found on RunFile

if (abs(IRC) == 1) then

  ! If this is the first macro iteration in the IRC search then
  ! pick up the reaction vector.

  if (Explicit_IRC .and. (iMEP == 0)) then
    ! Case 1)
    call dcopy_(3*nsAtom,TmpRx,1,MF,1)
  else if (iMEP == 0) then
    call NameRun('RUNOLD')
    call qpg_dArray('Reaction Vector',Found,nRx)
    !write(6,*) 'RUNOLD: Found=',Found
    if (Found) then
      ! Case 2)
      call Get_dArray('Reaction Vector',MF,3*nsAtom)
      call NameRun('#Pop')
    else
      call NameRun('#Pop')
      call qpg_dArray('Reaction Vector',Found,nRx)
      !write(6,*) 'RUNFILE: Found=',Found
      if (Found) then
        ! Case 3)
        call Get_dArray('Reaction Vector',MF,3*nsAtom)
      else
        call WarningMessage(2,'Error in RdCtl_Slapaf')
        write(6,*)
        write(6,*) '************ ERROR **************'
        write(6,*) 'IRC calculation but no IRC vector'
        write(6,*) '*********************************'
        call Quit_OnUserError()
      end if
    end if
  end if

  ! Fix the direction forward/backwards

  if ((iMEP == 0) .and. (iRC == -1)) call DScal_(3*nsAtom,-1.0d0,MF,1)
  if ((iMEP == 0) .and. (MEP_Type == 'TRANSVERSE')) call Put_dArray('Transverse',MF,3*nsAtom)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (FindTS .and. (.not. TSConstraints)) then
  call SysWarnMsg('RdCtl_Slapaf','WARNING:','FindTS specified, but no TSConstraints. '// &
                  'It is highly recommended to use TSConstraints in SLAPAF instead of (or in addition to) global constraints '// &
                  'when using FindTS. TSConstraints will be lifted in the final TS search.')
end if
TSConstraints = TSConstraints .and. FindTS

if ((MEP .or. rMEP) .and. (.not. Request_Alaska)) then

  ! If no initial direction given, use the gradient (force)

  call qpg_dArray('Transverse',Found,nRP)
  if (.not. Found) then
    ! Assume the initial reaction vector is already massaged
    if (Explicit_IRC) then
      call Put_dArray('Transverse',TmpRx,3*nsAtom)
    else
      ! The direction is given by the gradient, but in weighted coordinates
      call mma_allocate(Dir,3,nsAtom,Label='Dir')
      do iAtom=1,nsAtom
        xWeight = Weights(iAtom)
        do ixyz=1,3
          Dir(ixyz,iAtom) = Gx(ixyz,iAtom,iter)/xWeight
        end do
      end do
      call Put_dArray('Transverse',Dir,3*nsAtom)
      call mma_deallocate(Dir)
    end if
  end if
end if

if (Explicit_IRC) call mma_deallocate(TmpRx)

! Activate MPS update of Hessian if FindTS

if (FindTS) then

  if (iand(iOptH,64) == 64) then
    iOptH = ior(64,iand(iOptH,32)) ! EU
  else if (iand(iOptH,128) == 128) then
    iOptH = ior(128,iand(iOptH,32)) ! TS-BFGS
  else
    iOptH = ior(16,iand(iOptH,32)) ! MSP
  end if
  iOptC = ior(iOptC,4096)

  ! Increase the update window so that we will not lose the update
  ! which generated the negative curvature.

  nWndw = 4*nWndw
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify some options if TS search

if (iand(iOptC,128) /= 128) then
  if (iand(iOptH,8) /= 8) iOptH = ior(16,iand(iOptH,32)) ! MSP
  Line_search = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For TS optimization with the Saddle method set update to MSP

call qpg_dArray('Saddle',Found,nSaddle)
if (Found .and. (nSaddle /= 0)) then
  call mma_allocate(Tmp,nSaddle,Label='Tmp')
  call Get_dArray('Saddle',Tmp,nSaddle)
  HSR0 = Tmp(nSaddle-2)
  HSR = Tmp(nSaddle-1)
  Update = Tmp(nSaddle)
  if (Update == 2.0d0) then

    ! Enable FindTS procedure

    !write(6,*) 'Enable FindTS procedure'
    if (iand(iOptH,8) /= 8) iOptH = ior(16,iand(iOptH,32)) ! MSP
    nWndw = 4*nWndw
    ! make it look as if this were FindTS with constraints
    FindTS = .true.
    TSConstraints = .true.
    iOptC = ior(iOptC,4096)
    iOptC = ior(iOptC,8192)
    Beta = 0.1d0

  else

    ! Normal constrained optimization with a reduced threshold.
    ! Let the threshold be somewhat tighter as we are close to the TS.

    if ((HSR/HSR0 < 0.20d0) .or. (HSR < 0.20d0)) then
      !ThrGrd = 0.0003D0
      Beta = 0.1d0
    else
      !ThrGrd = 0.003D0
      ThrGrd = Ten*ThrGrd
    end if

    ! Add the constraints from the Saddle method

    call Merge_Constraints('UDC','UDC.Saddle','UDC',nLambda,iRow_c)

  end if
  call mma_deallocate(Tmp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nLbl = max(3*nsAtom+nLambda,iRow+iRow_c)
call mma_allocate(Lbl,nLbl,Label='Lbl')
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify some options if constraints are part of the calculation.

if ((nLambda > 0) .or. TSConstraints) then
  iOptC = ior(iOptC,256) ! Constraints
  Line_search = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! No iterations set iOptC=0

if (MxItr == 0) iOptC = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Activate some additional printing for numerical Hessian

!GGd: Coherency with patch 7.1.615 !      If (lNmHss) nPrint(122) = 10
!                                                                      *
!***********************************************************************
!                                                                      *
! Do some preprocessing due to input choice

if (Request_Alaska) nPrint(51) = 0
call PrePro(nsAtom,Cx(1,1,iter))
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of Kriging we use a sorting step in update_sl. For this
! to work we need the values of the internal coordinates for more
! points than the window size. Here we increase it with a factor of
! 2 temporarily. The sorted list will still be of the original size.
! However, the default window for kriging is twice as large as
! for conventional calculations.

if (Kriging) then
  nWndw = 4*nWndw  ! 2*2=4

  ! No micro iterations the first MEP iteration

  if ((MEP .or. rMEP) .and. (iter == 1)) Max_Microiterations = 0

  ! Reduce default maximum dispersion during the initial
  ! stage of a FindTS calculation: we don't want to fulfil the
  ! constraints too early

  call Qpg_iScalar('TS Search',Found)
  if (Found) call Get_lScalar('TS Search',Found)
  if (FindTS .and. (.not. (Found .or. Manual_Beta))) then
    Beta_Disp = 0.1d0
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out input parameters, No output if we didn't find the proper
! gradients on the runfile. We will come back!!!

if (SuperName == 'slapaf') then
  if (.not. Request_Alaska) call WrInp_sl()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
User_Def = iRow /= 0
if (.not. User_Def) nBVec = 1
!                                                                      *
!***********************************************************************
!                                                                      *
if (lNmHss .and. allocated(mRowH)) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '**************************************************'
  write(Lu,*) ' ERROR: NUMErical and ROWH are mutually exclusive '
  write(Lu,*) '**************************************************'
  call Quit_OnUserError()
end if
if (lCtoF .and. User_Def) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '******************************************'
  write(Lu,*) ' ERROR: CtoF and User-defined Coordinates '
  write(Lu,*) '        are mutually exclusive.           '
  write(Lu,*) '******************************************'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine RdCtl_Slapaf
