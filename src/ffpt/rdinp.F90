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
integer(kind=iwp) :: i, i1, i2, iCom, iEnd, iSta, j, k, newline, ctrl
real(kind=wp) :: C, W, X, Y, Z
integer(kind=iwp), parameter :: mCom = 11
character(len=*), parameter :: Com(mCom) = ['TITL','DIPO','EFLD','QUAD','OCTU','EFGR','RELA','GLBL','SELE','CUMU','END ']
logical(kind=iwp) :: Op0(9) = .false., Op2(3) = .false., Op3(4) = .false., Op4(8) = .false., Op5(12) = .false., Op6(8) = .false., &
                     skip

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

skip = .false.
mainLoop: do
  if (.not. skip) then
    Temp2 = Temp1
    Temp1 = Line
    read(u5,'(A)',iostat=ctrl) Line
    if (ctrl /= 0) call error(991)
    newline = newline+1
    Line = adjustl(Line)
    if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle mainLoop
    call StdFmt(Line,Token)
  end if
  skip = .true.

  select case (Token)
    case (Com(1))
      !---  Process the "TITL" command --------------------------------------*

      if (Op0(1)) call error(993)
      Op0(1) = .true.
      ComStk(1,0,0,0) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        mTit = mTit+1
        if (mTit <= MxTitL) then
          Title(mTit) = Line
        end if
      end do

    case (Com(2))
      !---  Process the "DIPO" command --------------------------------------*

      if (Op0(2)) call error(993)
      Op0(2) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        do i=1,len(Line)
          if (Line(i:i) == '=') Line(i:i) = ' '
          if (Line(i:i) == char(9)) Line(i:i) = ' '
        end do
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 2
        i2 = len(Line)
        select case (Token)
          case ('X')
            if (Op2(1)) call error(993)
            Op2(1) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,1,1,0) = .true.
            ComStk(2,1,1,1) = .true.
            ComVal(2,1,1,1) = W
          case ('Y')
            if (Op2(2)) call error(993)
            Op2(2) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,1,1,0) = .true.
            ComStk(2,1,1,2) = .true.
            ComVal(2,1,1,2) = W
          case ('Z')
            if (Op2(3)) call error(993)
            Op2(3) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,1,1,0) = .true.
            ComStk(2,1,1,3) = .true.
            ComVal(2,1,1,3) = W
          case default
            call error(992)
        end select
      end do

    case (Com(3))
      !---  Process the "EFLD" command --------------------------------------*

      if (Op0(3)) call error(993)
      Op0(3) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        do i=1,len(Line)
          if (Line(i:i) == '=') Line(i:i) = ' '
          if (Line(i:i) == char(9)) Line(i:i) = ' '
        end do
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 2
        i2 = len(Line)
        select case (Token)
          case ('X')
            if (Op3(1)) call error(993)
            Op3(1) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,3,1,0) = .true.
            ComStk(2,3,1,1) = .true.
            ComVal(2,3,1,1) = W
          case ('Y')
            if (Op3(2)) call error(993)
            Op3(2) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,3,1,0) = .true.
            ComStk(2,3,1,2) = .true.
            ComVal(2,3,1,2) = W
          case ('Z')
            if (Op3(3)) call error(993)
            Op3(3) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,3,1,0) = .true.
            ComStk(2,3,1,3) = .true.
            ComVal(2,3,1,3) = W
          case ('ORIG')
            if (Op3(4)) call error(993)
            Op3(4) = .true.
            i1 = 5
            read(Line(i1:i2),*,iostat=ctrl) X,Y,Z
            if (ctrl /= 0) call error(991)
            ComStk(2,3,2,0) = .true.
            ComStk(2,3,2,1) = .true.
            ComVal(2,3,2,1) = X
            ComStk(2,3,2,2) = .true.
            ComVal(2,3,2,2) = Y
            ComStk(2,3,2,3) = .true.
            ComVal(2,3,2,3) = Z
          case default
            call error(992)
        end select
      end do

    case (Com(4))
      !---  Process the "QUAD" command --------------------------------------*

      if (Op0(4)) call error(993)
      Op0(4) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (is_iostat_end(ctrl)) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        do i=1,len(Line)
          if (Line(i:i) == '=') Line(i:i) = ' '
          if (Line(i:i) == char(9)) Line(i:i) = ' '
        end do
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 3
        i2 = len(Line)
        select case (Token)
          case ('XX')
            if (Op4(1)) call error(993)
            Op4(1) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,1) = .true.
            ComVal(2,2,1,1) = W
          case ('YY')
            if (Op4(2)) call error(993)
            Op4(2) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,4) = .true.
            ComVal(2,2,1,4) = W
          case ('ZZ')
            if (Op4(3)) call error(993)
            Op4(3) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,6) = .true.
            ComVal(2,2,1,6) = W
          case ('XY')
            if (Op4(4)) call error(993)
            Op4(4) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,2) = .true.
            ComVal(2,2,1,2) = W
          case ('XZ')
            if (Op4(5)) call error(993)
            Op4(5) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,3) = .true.
            ComVal(2,2,1,3) = W
          case ('YZ')
            if (Op4(6)) call error(993)
            Op4(6) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,5) = .true.
            ComVal(2,2,1,5) = W
          case ('RR')
            if (Op4(7)) call error(993)
            Op4(7) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,2,1,0) = .true.
            ComStk(2,2,1,7) = .true.
            ComVal(2,2,1,7) = W
          case ('ORIG')
            if (Op4(8)) call error(993)
            Op4(8) = .true.
            i1 = 5
            read(Line(i1:i2),*,iostat=ctrl) X,Y,Z
            if (ctrl /= 0) call error(991)
            ComStk(2,2,2,0) = .true.
            ComStk(2,2,2,1) = .true.
            ComVal(2,2,2,1) = X
            ComStk(2,2,2,2) = .true.
            ComVal(2,2,2,2) = Y
            ComStk(2,2,2,3) = .true.
            ComVal(2,2,2,3) = Z
          case default
            call error(992)
        end select
      end do

    case (Com(5))
      !---  Process the "OCTU" command --------------------------------------*

      if (Op0(5)) call error(993)
      Op0(5) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        do i=1,len(Line)
          if (Line(i:i) == '=') Line(i:i) = ' '
          if (Line(i:i) == char(9)) Line(i:i) = ' '
        end do
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 4
        i2 = len(Line)
        select case (Token)
          case ('XXX')
            if (Op5(1)) call error(993)
            Op5(1) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,1) = .true.
            ComStk(2,6,1,1) = .true.
            ComVal(2,6,1,1) = W
          case ('XYY')
            if (Op5(2)) call error(993)
            Op5(2) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,4) = .true.
            ComStk(2,6,1,4) = .true.
            ComVal(2,6,1,4) = W
          case ('XZZ')
            if (Op5(3)) call error(993)
            Op5(3) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,6) = .true.
            ComStk(2,6,1,6) = .true.
            ComVal(2,6,1,6) = W
          case ('XXY')
            if (Op5(4)) call error(993)
            Op5(4) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,2) = .true.
            ComStk(2,6,1,2) = .true.
            ComVal(2,6,1,2) = W
          case ('YYY')
            if (Op5(5)) call error(993)
            Op5(5) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,7) = .true.
            ComStk(2,6,1,7) = .true.
            ComVal(2,6,1,7) = W
          case ('YZZ')
            if (Op5(6)) call error(993)
            Op5(6) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,9) = .true.
            ComStk(2,6,1,9) = .true.
            ComVal(2,6,1,9) = W
          case ('XXZ')
            if (Op5(7)) call error(993)
            Op5(8) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,3) = .true.
            ComStk(2,6,1,3) = .true.
            ComVal(2,6,1,3) = W
          case ('YYZ')
            if (Op5(9)) call error(993)
            Op5(9) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,8) = .true.
            ComStk(2,6,1,8) = .true.
            ComVal(2,6,1,8) = W
          case ('ZZZ')
            if (Op5(10)) call error(993)
            Op5(10) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,10) = .true.
            ComStk(2,6,1,10) = .true.
            ComVal(2,6,1,10) = W
          case ('XYZ')
            if (Op5(11)) call error(993)
            Op5(11) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,6,1,5) = .true.
            ComStk(2,6,1,5) = .true.
            ComVal(2,6,1,5) = W
          case ('ORIG')
            if (Op5(12)) call error(993)
            Op5(12) = .true.
            i1 = 5
            read(Line(i1:i2),*,iostat=ctrl) X,Y,Z
            if (ctrl /= 0) call error(991)
            ComStk(2,6,2,0) = .true.
            ComStk(2,6,2,1) = .true.
            ComVal(2,6,2,1) = X
            ComStk(2,6,2,2) = .true.
            ComVal(2,6,2,2) = Y
            ComStk(2,6,2,3) = .true.
            ComVal(2,6,2,3) = Z
          case default
            call error(992)
        end select
      end do

    case (Com(6))
      !---  Process the "EFGR" command --------------------------------------*

      if (Op0(6)) call error(993)
      Op0(6) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        do i=1,len(Line)
          if (Line(i:i) == '=') Line(i:i) = ' '
          if (Line(i:i) == char(9)) Line(i:i) = ' '
        end do
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 3
        i2 = len(Line)
        select case (Token)
          case ('XX')
            if (Op6(1)) call error(993)
            Op6(1) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,1) = .true.
            ComVal(2,4,1,1) = W
          case ('YY')
            if (Op6(2)) call error(993)
            Op6(2) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,4) = .true.
            ComVal(2,4,1,4) = W
          case ('ZZ')
            if (Op6(3)) call error(993)
            Op6(3) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,6) = .true.
            ComVal(2,4,1,6) = W
          case ('XY')
            if (Op6(4)) call error(993)
            Op6(4) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,2) = .true.
            ComVal(2,4,1,2) = W
          case ('XZ')
            if (Op6(5)) call error(993)
            Op6(5) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,3) = .true.
            ComVal(2,4,1,3) = W
          case ('YZ')
            if (Op6(6)) call error(993)
            Op6(6) = .true.
            read(Line(i1:i2),*,iostat=ctrl) W
            if (ctrl /= 0) call error(991)
            ComStk(2,4,1,0) = .true.
            ComStk(2,4,1,5) = .true.
            ComVal(2,4,1,5) = W
          case ('ORIG')
            if (Op6(8)) call error(993)
            Op6(8) = .true.
            i1 = 5
            read(Line(i1:i2),*,iostat=ctrl) X,Y,Z
            if (ctrl /= 0) call error(991)
            ComStk(2,4,2,0) = .true.
            ComStk(2,4,2,1) = .true.
            ComVal(2,4,2,1) = X
            ComStk(2,4,2,2) = .true.
            ComVal(2,4,2,2) = Y
            ComStk(2,4,2,3) = .true.
            ComVal(2,4,2,3) = Z
          case default
            call error(992)
        end select
      end do

    case (Com(7))
      !---  Process the "RELA" command --------------------------------------*

      if (Op0(7)) call error(993)
      Op0(7) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (.not. (Line(1:1) == ' ' .or. Line(1:1) == '*')) exit
      end do
      i1 = 1
      i2 = len(Line)
      read(Line(i1:i2),*,iostat=ctrl) W
      if (ctrl /= 0) call error(991)
      ComStk(2,5,0,0) = .true.
      ComStk(2,5,0,1) = .true.
      ComVal(2,5,0,1) = W
      skip = .false.
      cycle mainLoop

    case (Com(8))
      !---  Process the "GLBL" command --------------------------------------*

      if (Op0(8)) call error(993)
      Op0(8) = .true.
      ComStk(3,0,0,0) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (Line(1:1) == ' ' .or. Line(1:1) == '*') cycle
        call UpCase(Line)
        call StdFmt(Line,Token)
        do iCom=1,mCom
          if (Token == Com(iCom)) cycle mainLoop
        end do
        i1 = 1
        i2 = len(Line)
        mLbl = mLbl+1
        if (mLbl > MxLbl) call error(994)
        iSta = index(Line,'''')
        Line(i1:iSta) = ' '
        iEnd = index(Line,'''')
        gLblN(mLbl) = Line(iSta+1:iEnd-1)
        Line(i1:iEnd) = ' '
        read(Line(i1:i2),*,iostat=ctrl) C,W
        if (ctrl /= 0) call error(991)
        gLblC(mLbl) = nint(C,kind=iwp)
        gLblW(mLbl) = W
        do i=1,mLbl-1
          if (gLblN(i) == gLblN(mLbl) .and. gLblC(i) == gLblC(mLbl)) call error(993)
        end do
      end do

    case (Com(9))
      !---  Process the "SELE" command --------------------------------------*

      !-- Initialize
      ComStk(4,0,0,0) = .true.
      do
        Temp2 = Temp1
        Temp1 = Line
        newline = newline+1
        read(u5,'(A)',iostat=ctrl) Line
        if (ctrl /= 0) call error(991)
        Line = adjustl(Line)
        if (.not. (Line(1:1) == ' ' .or. Line(1:1) == '*')) exit
      end do
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
        if (i < 10) write(FmtLog,'(A,I1,A)') '(',i-1,'L2)'
        if (i >= 10) write(FmtLog,'(A,I2,A)') '(',i-1,'L2)'
        read(u5,FmtLog) (Bonds(i,j),j=1,i-1)
        do j=1,i-1
          Bonds(j,i) = Bonds(i,j)
        end do
      end do
      read(u5,*) (TranCoo(k),k=1,3)
      !read(u5,*) SiffBond
      skip = .false.
      cycle mainLoop

    case (Com(10))
      !---  Process the "CUMU" command --------------------------------------*

      ! Add the perturbation to the current H0 instead of
      ! to the vacuum H0. This enables multiple FFPT runs after eachother,
      ! useful when using SELE.
      LCumulate = .true.
      skip = .false.
      cycle mainLoop

    case (Com(11))
      !---  Process the "END " command --------------------------------------*

      if (Op0(9)) call error(993)
      Op0(9) = .true.
      !mStk(4,0,0,0) = .true.
      ComStk(5,0,0,0) = .true.
      !---  Check for redunancy in the origin input
      if (ComStk(2,2,2,4) .and. (ComStk(2,2,2,1) .or. ComStk(2,2,2,2) .or. ComStk(2,2,2,3))) call error(996)
      if (ComStk(2,3,2,4) .and. (ComStk(2,3,2,1) .or. ComStk(2,3,2,2) .or. ComStk(2,3,2,3))) call error(996)
      if (ComStk(2,4,2,4) .and. (ComStk(2,4,2,1) .or. ComStk(2,4,2,2) .or. ComStk(2,4,2,3))) call error(996)
      if (ComStk(2,6,2,4) .and. (ComStk(2,6,2,1) .or. ComStk(2,6,2,2) .or. ComStk(2,6,2,3))) call error(996)
      exit mainLoop

    case default
      call error(992)
  end select
end do mainLoop

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

contains

!----------------------------------------------------------------------*
!     Error handling                                                   *
!----------------------------------------------------------------------*
subroutine error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (991)
      write(u6,*)
      write(u6,'(2X,A)') 'The program failed to read the input.'
      write(u6,'(2X,A)') 'Please check your input data.'
      write(u6,*)
      write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
      write(u6,'(2X,A,A)') 'The current line is:      ',Line
      write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
      write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
      call Quit_OnUserError()
    case (992)
      write(u6,*)
      write(u6,'(2X,A)') 'The program has been supplied with an unknown'
      write(u6,'(2X,A)') 'keyword. Please correct your input data.'
      write(u6,*)
      write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
      write(u6,'(2X,A,A)') 'The current line is:      ',Line
      write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
      write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
      call Quit_OnUserError()
    case (993)
      write(u6,*)
      write(u6,'(2X,A)') 'A command or one of its components has been'
      write(u6,'(2X,A)') 'multiply defined. Please correct your input.'
      write(u6,*)
      write(u6,'(2X,A,I3.3,A)') 'The error occured at line',newline,' after the &FFPT &END line'
      write(u6,'(2X,A,A)') 'The current line is:      ',Line
      write(u6,'(2X,A,A)') 'The previous line is:     ',Temp1
      write(u6,'(2X,A,A)') 'The next previous line is:',Temp2
      call Quit_OnUserError()
    case (994)
      write(u6,*)
      write(u6,'(2X,A)') 'The number of perturbations requested exceeds'
      write(u6,'(2X,A)') 'the internal buffer size. Increase the para-'
      write(u6,'(2X,A)') 'meter MxLbl and recompile the program.'
      call Quit_OnUserError()
    case (996)
      write(u6,*)
      write(u6,'(2X,A)') 'The definition of the origin of an operator '
      write(u6,'(2X,A)') 'is not unique.  Please correct your input.'
    case default
  end select

  call Quit_OnUserError()

end subroutine error

end subroutine RdInp_FFPT
