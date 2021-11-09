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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!               1999, Niclas Forsberg                                  *
!               Per Ake Malmqvist                                      *
!***********************************************************************

subroutine MulaRdNLst(iUnit,NameIn)
!  Locate the beginning of an input stream
!  (similar to FORTRAN NAMELIST read known to some systems)
!
!  Calling arguments:
!    iUnit  : Type integer, input
!             FORTRAN unit number
!    NameIn : Type character string, input
!             Character string marking the beginning of the input
!
!  Written by:
!    M.P. Fuelscher and P.O. Widmark
!    University of Lund, Sweden, 1993
!
!  Slightly modified by:
!    Niclas Forsberg, 1999

use Definitions, only: u6

#include "inputdata.fh"
character NameIn*(*)
character*8 StdNam
character*80 Line
integer StrnLn

! Convert the Name to internal standard format.
call StdFmt(NameIn,StdNam)
lStdNam = StrnLn(StdNam)

! Read until an input Line is located which starts with
! the string, Name, not before the second column
do while (.true.)
  read(iUnit,'(A80)',end=900) Line
  call LeftAd(Line)
  call UpCase(Line)
  if ((Line(1:1) == '&') .and. (Line(2:lStdNam+1) == StdNam(1:lStdNam))) then
    return
  end if
end do

! Error exit
900 continue
write(u6,*) 'MulaRdNLst error: Could not locate input.'
write(u6,*) 'Looking for: &'//StdNam
call Quit_OnUserError()

end subroutine MulaRdNLst
!####
subroutine WordPos(k,InLine,iStart,iStop)
!  Purpose:
!    Return the position (iStart,iStop) of the first word after
!    position k in InLine.
!
!  Input:
!    k        : Integer - the first position after keyword in the
!               string InLine.
!    InLine   : Character string.
!
!  Output:
!    iStart   : Integer - the position of the first non-blank
!               character after position k in InLine.
!    iStop,k  : Integer - the position of the first blank character
!               after position iStart in InLine.

#include "inputdata.fh"
character*1 Ch
character InLine*(*)

Length = len(InLine)
Ch = InLine(k:k)
do while ((Ch == ' ') .and. (k < Length))
  k = k+1
  Ch = InLine(k:k)
end do
if (k+1 >= Length) then
  iStart = Length
  iStop = Length
  return
end if
iStart = k
k = k+1
Ch = InLine(k:k)
do while ((Ch /= ' ') .and. (k < Length))
  k = k+1
  Ch = InLine(k:k)
end do
iStop = k-1

end subroutine WordPos
!####
subroutine KeyWord(nUnit,KeyWd,rewind,Exist)
!  Purpose:
!    Read from a file until a line contains the string KeyWd.
!
!  Input:
!    nUnit     : Unit to read from.
!    KeyWd     : Character string - keyword.
!    Rewind    : Logical.
!
!  Output:
!    File      : Next line in the file contains the data connected
!                with the keyword KeyWd.
!  Calls:
!    MulaRdNlst
!    Normalize

use Definitions, only: u6

#include "inputdata.fh"
character KeyWd*(*)
character*32 TmpKey
character*80 InLine, OutLine
integer i, KLen, BlInd
logical rewind, Exist

InLine = KeyWd
call Normalize(InLine,OutLine)
! Index of first blank character from the end:
BlInd = 81
do i=80,1,-1
  if (Outline(i:i) /= ' ') goto 15
  BlInd = i
end do
15 continue
if (BlInd > 80) then
  write(u6,*) 'KEYWORD: KeyWd is too long.'
  write(u6,*) 'KeyWd:',KeyWd
  write(u6,*) 'After normalization (OutLine):'
  write(u6,*) OutLine
  call abend()
end if
Exist = .false.
if (BlInd == 1) return

KLen = BlInd-1
TmpKey = OutLine(1:KLen)

if (rewind) then
  rewind nUnit
  call MulaRdNlst(nUnit,InputName)
end if

read(nUnit,'(a80)',end=10,err=20) InLine
call Normalize(InLine,OutLine)
do while (OutLine(1:KLen) /= TmpKey(1:KLen))
  read(nUnit,'(a80)',end=10,err=20) InLine
  call Normalize(InLine,OutLine)
end do
10 continue
Exist = OutLine(1:KLen) == TmpKey(1:KLen)

return

20 write(u6,*) ' I/O error on unit nUnit=',nUnit
call abend()

end subroutine KeyWord
!####
real*8 function StrToDble(InString)
!  Purpose:
!    Convert a number in string InString to Real*8.
!
!  Input:
!    InString   : Character string.
!
!  Output:
!    xNum     : Real*8.

use Definitions, only: u6

#include "inputdata.fh"
real*8 xNum
integer length
character InString*(*)
character*1 OneDigit
character*2 TwoDigits

length = len(InString)
write(u6,*) ' The string is:',InString
write(u6,*) ' Its length is:',length
if (length < 10) then
  write(OneDigit,'(i1)') length
  write(u6,*) '   OneDigit is:',OneDigit
  write(u6,*) ' The format is:','(f'//OneDigit//'.0)'
  read(InString,'(f'//OneDigit//'.0)') xnum
else
  write(TwoDigits,'(i2)') length
  write(u6,*) '  TwoDigits is:',TwoDigits
  write(u6,*) ' The format is:','(f'//TwoDigits//'.0)'
  read(InString,'(f'//TwoDigits//'.0)') xNum
end if
StrToDble = xNum

end function StrToDble
!####
function iStrToInt(InLine)
!  Purpose:
!    Convert a number in string InLine to integer.
!
!  Input:
!    InLine   : Character string.
!
!  Output:
!    num      : Integer.

integer num
character InLine*(*)
character ch
logical minus

IntVal(ch) = index('0123456789',ch)-1

length = len(InLine)

! Convert from string to integer.
isum = 0
i = 0
minus = .false.
do nPos=length,1,-1
  if (InLine(nPos:nPos) == '-') then
    minus = .true.
  else
    ch = InLine(nPos:nPos)
    isum = isum+IntVal(ch)*10**(i)
    i = i+1
  end if
end do

if (minus) then
  num = -isum
else
  num = isum
end if

iStrToInt = num

end function iStrToInt
!####
subroutine Normalize(line,line2)
!  Purpose:
!    Convert a string to uppercase characters.
!
!  Input:
!    line     : Character string.
!
!  Output:
!    line2    : Character string - uppercase letters.
!
!  Written by:
!    P-AA Malmquist,
!    Dept. of Theoretical Chemistry, Lund University.

use Definitions, only: u6

parameter(ntrans=256)
#include "inputdata.fh"
character line*(*)
character line2*(*)
character*26 lcase, ucase
character*1 ch
integer itrans(ntrans)
logical wrdend
intrinsic ichar, char, len
data icalled/0/
data lcase,ucase/'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

! Set up character transformation table:
if (icalled == 0) then
  do i=1,ntrans
    itrans(i) = i
  end do
  do j=1,len(lcase)
    i = ichar(lcase(j:j))
    itrans(i) = ichar(ucase(j:j))
  end do
end if
len1 = len(line)
len2max = len(line2)
ista = 0
! Find first non-blank character, at position ista.
do i=1,len1
  ch = line(i:i)
  if (ch /= ' ') then
    ista = i
    goto 11
  end if
end do
11 continue
line2 = ' '
if (ista == 0) return
wrdend = .false.
len2 = 0
do i=ista,len1
  ch = line(i:i)
  if (ch /= ' ') then
    if (wrdend) then
      if (len2 > len2max-1) goto 99
      line2(len2+1:len2+1) = ' '
      len2 = len2+1
    end if
    if (len2 > len2max-1) goto 99
    ch = char(itrans(ichar(ch)))
    line2(len2+1:len2+1) = ch
    len2 = len2+1
    wrdend = .false.
  else
    wrdend = .true.
  end if
end do

return

99 continue
write(u6,*) ' WARNING: Line truncated in NORMALIZE.'

end subroutine Normalize
