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
! Copyright (C) 2000, Markus P. Fuelscher                              *
!***********************************************************************

subroutine RdInp_FFPT()
!***********************************************************************
!                                                                      *
!     Objective: Read and interprete input                             *
!                                                                      *
!***********************************************************************
!                                                                      *
!     In order to introduce the ffpt input to the CERIUS2 interface    *
!     the parser had to be changed!                                    *
!     M. Fuelscher, Lund Univeristy, Sweden, February 2000             *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: Atoms, TranCoo, LCumulate, iSelection, Bonds, mLbl, mTit, MxTitL, nSets, MxLbl, Title, ComStk, ComVal, &
                       gLblN, gLblC, gLblW
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u5, u6

implicit none
character(len=20) :: FmtLog
character(len=72) :: Line, Temp1, Temp2
character(len=4) :: Token
integer(kind=iwp) :: i, i1, i2, iCom, iEnd, iSta, j, jCom, k, newline
real(kind=wp) :: C, W, X, Y, Z
integer(kind=iwp), parameter :: mCom = 11
character(len=4), parameter :: Com(mCom) = ['TITL','DIPO','EFLD','QUAD','OCTU','EFGR','RELA','GLBL','SELE','CUMU','END ']
logical(kind=iwp) :: Op0(9) = .false., Op2(3) = .false., Op3(4) = .false., Op4(8) = .false., Op5(12) = .false., Op6(8) = .false.

LCumulate = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Locate "start of input"                                          *
!                                                                      *
!----------------------------------------------------------------------*

call RdNlst(u5,'FFPT')
Temp2 = ' '
Temp1 = ' '
Line = ' &FFPT &END'

!----------------------------------------------------------------------*
!     Initialize counters                                              *
!----------------------------------------------------------------------*

newline = 0
mLbl = 0
mTit = 0

!----------------------------------------------------------------------*
!     Read the input stream line by line and identify key command      *
!----------------------------------------------------------------------*

jCom = 0
1 Temp2 = Temp1
Temp1 = Line
read(u5,'(A)',Err=991,end=991) Line
newline = newline+1
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 1
call StdFmt(Line,Token)
2 jCom = 0
do iCom=1,mCom
  if (Token == Com(iCom)) jCom = iCom
end do
if (jCom == 0) goto 992

!----------------------------------------------------------------------*
!     Branch to the processing of the command sections                 *
!----------------------------------------------------------------------*

goto(10,20,30,40,50,60,70,80,90,100,1000),jCom

!---  Process the "TITL" command --------------------------------------*

10 continue
if (Op0(1)) goto 993
Op0(1) = .true.
ComStk(1,0,0,0) = .true.
15 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 15
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
mTit = mTit+1
if (mTit <= MxTitL) then
  Title(mTit) = Line
end if
goto 15

!---  Process the "DIPO" command --------------------------------------*

20 continue
if (Op0(2)) goto 993
Op0(2) = .true.
25 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 25
call UpCase(Line)
do i=1,len(Line)
  if (Line(i:i) == '=') Line(i:i) = ' '
  if (Line(i:i) == char(9)) Line(i:i) = ' '
end do
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 2
i2 = len(Line)
if (Token == 'X') then
  if (Op2(1)) goto 993
  Op2(1) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,1,1,0) = .true.
  ComStk(2,1,1,1) = .true.
  ComVal(2,1,1,1) = W
  goto 25
end if
if (Token == 'Y') then
  if (Op2(2)) goto 993
  Op2(2) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,1,1,0) = .true.
  ComStk(2,1,1,2) = .true.
  ComVal(2,1,1,2) = W
  goto 25
end if
if (Token == 'Z') then
  if (Op2(3)) goto 993
  Op2(3) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,1,1,0) = .true.
  ComStk(2,1,1,3) = .true.
  ComVal(2,1,1,3) = W
  goto 25
end if
goto 992

!---  Process the "EFLD" command --------------------------------------*

30 continue
if (Op0(3)) goto 993
Op0(3) = .true.
35 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 35
call UpCase(Line)
do i=1,len(Line)
  if (Line(i:i) == '=') Line(i:i) = ' '
  if (Line(i:i) == char(9)) Line(i:i) = ' '
end do
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 2
i2 = len(Line)
if (Token == 'X') then
  if (Op3(1)) goto 993
  Op3(1) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,3,1,0) = .true.
  ComStk(2,3,1,1) = .true.
  ComVal(2,3,1,1) = W
  goto 35
end if
if (Token == 'Y') then
  if (Op3(2)) goto 993
  Op3(2) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,3,1,0) = .true.
  ComStk(2,3,1,2) = .true.
  ComVal(2,3,1,2) = W
  goto 35
end if
if (Token == 'Z') then
  if (Op3(3)) goto 993
  Op3(3) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,3,1,0) = .true.
  ComStk(2,3,1,3) = .true.
  ComVal(2,3,1,3) = W
  goto 35
end if
i1 = 5
if (Token == 'ORIG') then
  if (Op3(4)) goto 993
  Op3(4) = .true.
  i1 = 5
  read(Line(i1:i2),*,Err=991,end=991) X,Y,Z
  ComStk(2,3,2,0) = .true.
  ComStk(2,3,2,1) = .true.
  ComVal(2,3,2,1) = X
  ComStk(2,3,2,2) = .true.
  ComVal(2,3,2,2) = Y
  ComStk(2,3,2,3) = .true.
  ComVal(2,3,2,3) = Z
  goto 35
end if
goto 992

!---  Process the "QUAD" command --------------------------------------*

40 continue
if (Op0(4)) goto 993
Op0(4) = .true.
45 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 45
call UpCase(Line)
do i=1,len(Line)
  if (Line(i:i) == '=') Line(i:i) = ' '
  if (Line(i:i) == char(9)) Line(i:i) = ' '
end do
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 3
i2 = len(Line)
if (Token == 'XX') then
  if (Op4(1)) goto 993
  Op4(1) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,1) = .true.
  ComVal(2,2,1,1) = W
  goto 45
end if
if (Token == 'YY') then
  if (Op4(2)) goto 993
  Op4(2) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,4) = .true.
  ComVal(2,2,1,4) = W
  goto 45
end if
if (Token == 'ZZ') then
  if (Op4(3)) goto 993
  Op4(3) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,6) = .true.
  ComVal(2,2,1,6) = W
  goto 45
end if
if (Token == 'XY') then
  if (Op4(4)) goto 993
  Op4(4) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,2) = .true.
  ComVal(2,2,1,2) = W
  goto 45
end if
if (Token == 'XZ') then
  if (Op4(5)) goto 993
  Op4(5) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,3) = .true.
  ComVal(2,2,1,3) = W
  goto 45
end if
if (Token == 'YZ') then
  if (Op4(6)) goto 993
  Op4(6) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,5) = .true.
  ComVal(2,2,1,5) = W
  goto 45
end if
if (Token == 'RR') then
  if (Op4(7)) goto 993
  Op4(7) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,2,1,0) = .true.
  ComStk(2,2,1,7) = .true.
  ComVal(2,2,1,7) = W
  goto 45
end if
if (Token == 'ORIG') then
  if (Op4(8)) goto 993
  Op4(8) = .true.
  i1 = 5
  read(Line(i1:i2),*,Err=991,end=991) X,Y,Z
  ComStk(2,2,2,0) = .true.
  ComStk(2,2,2,1) = .true.
  ComVal(2,2,2,1) = X
  ComStk(2,2,2,2) = .true.
  ComVal(2,2,2,2) = Y
  ComStk(2,2,2,3) = .true.
  ComVal(2,2,2,3) = Z
  goto 45
end if
goto 992

!---  Process the "OCTU" command --------------------------------------*

50 continue
if (Op0(5)) goto 993
Op0(5) = .true.
55 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 55
call UpCase(Line)
do i=1,len(Line)
  if (Line(i:i) == '=') Line(i:i) = ' '
  if (Line(i:i) == char(9)) Line(i:i) = ' '
end do
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 4
i2 = len(Line)
if (Token == 'XXX') then
  if (Op5(1)) goto 993
  Op5(1) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,1) = .true.
  ComStk(2,6,1,1) = .true.
  ComVal(2,6,1,1) = W
  goto 55
end if
if (Token == 'XYY') then
  if (Op5(2)) goto 993
  Op5(2) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,4) = .true.
  ComStk(2,6,1,4) = .true.
  ComVal(2,6,1,4) = W
  goto 55
end if
if (Token == 'XZZ') then
  if (Op5(3)) goto 993
  Op5(3) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,6) = .true.
  ComStk(2,6,1,6) = .true.
  ComVal(2,6,1,6) = W
  goto 55
end if
if (Token == 'XXY') then
  if (Op5(4)) goto 993
  Op5(4) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,2) = .true.
  ComStk(2,6,1,2) = .true.
  ComVal(2,6,1,2) = W
  goto 55
end if
if (Token == 'YYY') then
  if (Op5(5)) goto 993
  Op5(5) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,7) = .true.
  ComStk(2,6,1,7) = .true.
  ComVal(2,6,1,7) = W
  goto 55
end if
if (Token == 'YZZ') then
  if (Op5(6)) goto 993
  Op5(6) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,9) = .true.
  ComStk(2,6,1,9) = .true.
  ComVal(2,6,1,9) = W
  goto 55
end if
if (Token == 'XXZ') then
  if (Op5(7)) goto 993
  Op5(8) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,3) = .true.
  ComStk(2,6,1,3) = .true.
  ComVal(2,6,1,3) = W
  goto 55
end if
if (Token == 'YYZ') then
  if (Op5(9)) goto 993
  Op5(9) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,8) = .true.
  ComStk(2,6,1,8) = .true.
  ComVal(2,6,1,8) = W
  goto 55
end if
if (Token == 'ZZZ') then
  if (Op5(10)) goto 993
  Op5(10) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,10) = .true.
  ComStk(2,6,1,10) = .true.
  ComVal(2,6,1,10) = W
  goto 55
end if
if (Token == 'XYZ') then
  if (Op5(11)) goto 993
  Op5(11) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,6,1,5) = .true.
  ComStk(2,6,1,5) = .true.
  ComVal(2,6,1,5) = W
  goto 55
end if
if (Token == 'ORIG') then
  if (Op5(12)) goto 993
  Op5(12) = .true.
  i1 = 5
  read(Line(i1:i2),*,Err=991,end=991) X,Y,Z
  ComStk(2,6,2,0) = .true.
  ComStk(2,6,2,1) = .true.
  ComVal(2,6,2,1) = X
  ComStk(2,6,2,2) = .true.
  ComVal(2,6,2,2) = Y
  ComStk(2,6,2,3) = .true.
  ComVal(2,6,2,3) = Z
  goto 55
end if
goto 992

!---  Process the "EFGR" command --------------------------------------*

60 continue
if (Op0(6)) goto 993
Op0(6) = .true.
65 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 65
call UpCase(Line)
do i=1,len(Line)
  if (Line(i:i) == '=') Line(i:i) = ' '
  if (Line(i:i) == char(9)) Line(i:i) = ' '
end do
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 3
i2 = len(Line)
if (Token == 'XX') then
  if (Op6(1)) goto 993
  Op6(1) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,1) = .true.
  ComVal(2,4,1,1) = W
  goto 65
end if
if (Token == 'YY') then
  if (Op6(2)) goto 993
  Op6(2) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,4) = .true.
  ComVal(2,4,1,4) = W
  goto 65
end if
if (Token == 'ZZ') then
  if (Op6(3)) goto 993
  Op6(3) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,6) = .true.
  ComVal(2,4,1,6) = W
  goto 65
end if
if (Token == 'XY') then
  if (Op6(4)) goto 993
  Op6(4) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,2) = .true.
  ComVal(2,4,1,2) = W
  goto 65
end if
if (Token == 'XZ') then
  if (Op6(5)) goto 993
  Op6(5) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,3) = .true.
  ComVal(2,4,1,3) = W
  goto 65
end if
if (Token == 'YZ') then
  if (Op6(6)) goto 993
  Op6(6) = .true.
  read(Line(i1:i2),*,Err=991,end=991) W
  ComStk(2,4,1,0) = .true.
  ComStk(2,4,1,5) = .true.
  ComVal(2,4,1,5) = W
  goto 65
end if
if (Token == 'ORIG') then
  if (Op6(8)) goto 993
  Op6(8) = .true.
  i1 = 5
  read(Line(i1:i2),*,Err=991,end=991) X,Y,Z
  ComStk(2,4,2,0) = .true.
  ComStk(2,4,2,1) = .true.
  ComVal(2,4,2,1) = X
  ComStk(2,4,2,2) = .true.
  ComVal(2,4,2,2) = Y
  ComStk(2,4,2,3) = .true.
  ComVal(2,4,2,3) = Z
  goto 65
end if
goto 992

!---  Process the "RELA" command --------------------------------------*

70 continue
if (Op0(7)) goto 993
Op0(7) = .true.
75 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 75
i1 = 1
i2 = len(Line)
read(Line(i1:i2),*,Err=991,end=991) W
ComStk(2,5,0,0) = .true.
ComStk(2,5,0,1) = .true.
ComVal(2,5,0,1) = W
goto 1

!---  Process the "GLBL" command --------------------------------------*

80 continue
if (Op0(8)) goto 993
Op0(8) = .true.
ComStk(3,0,0,0) = .true.
85 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 85
call UpCase(Line)
call StdFmt(Line,Token)
do iCom=1,mCom
  if (Token == Com(iCom)) goto 2
end do
i1 = 1
i2 = len(Line)
mLbl = mLbl+1
if (mLbl > MxLbl) goto 994
iSta = index(Line,'''')
Line(i1:iSta) = ' '
iEnd = index(Line,'''')
gLblN(mLbl) = Line(iSta+1:iEnd-1)
Line(i1:iEnd) = ' '
read(Line(i1:i2),*,end=991) C,W
gLblC(mLbl) = nint(C)
gLblW(mLbl) = W
do i=1,mLbl-1
  if (gLblN(i) == gLblN(mLbl) .and. gLblC(i) == gLblC(mLbl)) goto 993
end do
goto 85

!---  Process the "SELE" command --------------------------------------*

90 continue
!-- Initialize
ComStk(4,0,0,0) = .true.
95 Temp2 = Temp1
Temp1 = Line
newline = newline+1
read(u5,'(A)',Err=991,end=991) Line
call LeftAd(Line)
if (Line(1:1) == ' ' .or. Line(1:1) == '*') goto 95
read(Line,*) nSets
call mma_allocate(iSelection,2,nSets,label='iSelection')
call mma_allocate(Atoms,nSets,label='Atoms')
call mma_allocate(Bonds,nSets,nSets,label='Bonds')
Atoms(:) = .false.
Bonds(:,:) = .false.
do i=1,nSets
  read(u5,*) Atoms(i),iSelection(1,i),iSelection(2,i)
end do
do i=2,nSets
  if (i < 10) write(FmtLog,79121) '(',i-1,'L2)'
  if (i >= 10) write(FmtLog,79122) '(',i-1,'L2)'
  read(u5,FmtLog) (Bonds(i,j),j=1,i-1)
  do j=1,i-1
    Bonds(j,i) = Bonds(i,j)
  end do
end do
read(u5,*) (TranCoo(k),k=1,3)
!read(u5,*) SiffBond
79121 format(A,I1,A)
79122 format(A,I2,A)
Go to 1

!---  Process the "CUMU" command --------------------------------------*

100 continue
!     Add the perturbation to the current H0 instead of
!     to the vacuum H0. This enables multiple FFPT runs after eachother,
!     useful when using SELE.
LCumulate = .true.
Go to 1

!---  Process the "END " command --------------------------------------*

1000 continue
if (Op0(9)) goto 993
Op0(9) = .true.
!mStk(4,0,0,0) = .true.
ComStk(5,0,0,0) = .true.
!---  Check for redunancy in the origin input
if (ComStk(2,2,2,4)) then
  if (ComStk(2,2,2,1) .or. ComStk(2,2,2,2) .or. ComStk(2,2,2,3)) goto 996
end if
if (ComStk(2,3,2,4)) then
  if (ComStk(2,3,2,1) .or. ComStk(2,3,2,2) .or. ComStk(2,3,2,3)) goto 996
end if
if (ComStk(2,4,2,4)) then
  if (ComStk(2,4,2,1) .or. ComStk(2,4,2,2) .or. ComStk(2,4,2,3)) goto 996
end if
if (ComStk(2,6,2,4)) then
  if (ComStk(2,6,2,1) .or. ComStk(2,6,2,2) .or. ComStk(2,6,2,3)) goto 996
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

!----------------------------------------------------------------------*
!     Error handling                                                   *
!----------------------------------------------------------------------*

991 write(u6,*)
write(u6,'(2X,A)') 'The program failed to read the input.'
write(u6,'(2X,A)') 'Please check your input data.'
write(u6,*)
write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
write(u6,'(2X,A,A)') 'The current line is:      ',Line
write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
call Quit_OnUserError()
992 write(u6,*)
write(u6,'(2X,A)') 'The program has been supplied with an unknown'
write(u6,'(2X,A)') 'keyword. Please correct your input data.'
write(u6,*)
write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
write(u6,'(2X,A,A)') 'The current line is:      ',Line
write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
call Quit_OnUserError()
993 write(u6,*)
write(u6,'(2X,A)') 'A command or one of its components has been'
write(u6,'(2X,A)') 'multiply defined. Please correct your input.'
write(u6,*)
write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
write(u6,'(2X,A,A)') 'The current line is:      ',Line
write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
call Quit_OnUserError()
994 write(u6,*)
write(u6,'(2X,A)') 'The number of perturbations requested exceeds'
write(u6,'(2X,A)') 'the internal buffer size. Increase the para-'
write(u6,'(2X,A)') 'meter MxLbl and recompile the program.'
call Quit_OnUserError()
996 write(u6,*)
write(u6,'(2X,A)') 'The definition of the origin of an operator '
write(u6,'(2X,A)') 'is not unique.  Please correct your input.'
call Quit_OnUserError()
end subroutine RdInp_FFPT
