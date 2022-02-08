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

!module IOTools

!contains

subroutine MulaRdNLst(iUnit,NameIn)
!  Locate the beginning of an input stream
!  (similar to FORTRAN NAMELIST read known to some systems)
!
!  Calling arguments:
!    iUnit  : Type integer, input - FORTRAN unit number
!    NameIn : Type character string, input - Character string marking the beginning of the input
!
!  Written by:
!    M.P. Fuelscher and P.O. Widmark
!    University of Lund, Sweden, 1993
!
!  Slightly modified by:
!    Niclas Forsberg, 1999

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iUnit
character(len=*), intent(in) :: NameIn
integer(kind=iwp) :: istatus, lStdNam
character(len=80) :: Line
character(len=8) :: StdNam
integer(kind=iwp), external :: StrnLn

! Convert the Name to internal standard format.
call StdFmt(NameIn,StdNam)
lStdNam = StrnLn(StdNam)

! Read until an input Line is located which starts with
! the string, Name, not before the second column
do while (.true.)
  read(iUnit,'(A80)',iostat=istatus) Line
  if (istatus < 0) then
    write(u6,*) 'MulaRdNLst error: Could not locate input.'
    write(u6,*) 'Looking for: &'//StdNam
    call Quit_OnUserError()
  end if
  call UpCase(Line)
  Line = adjustl(Line)
  if ((Line(1:1) == '&') .and. (Line(2:lStdNam+1) == StdNam(1:lStdNam))) exit
end do

end subroutine MulaRdNLst

subroutine WordPos(k,InLine,iStart,iStop)
!  Purpose:
!    Return the position (iStart,iStop) of the first word after position k in InLine.
!
!  Input:
!    k        : Integer - the first position after keyword in the string InLine.
!    InLine   : Character string.
!
!  Output:
!    iStart   : Integer - the position of the first non-blank character after position k in InLine.
!    iStop,k  : Integer - the position of the first blank character after position iStart in InLine.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: k
character(len=*), intent(in) :: InLine
integer(kind=iwp), intent(out) :: iStart, iStop
integer(kind=iwp) :: Length
character :: Ch

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

subroutine KeyWord(nUnit,KeyWd,rwnd,Exists)
!  Purpose:
!    Read from a file until a line contains the string KeyWd.
!
!  Input:
!    nUnit     : Unit to read from.
!    KeyWd     : Character string - keyword.
!    rwnd      : Logical.
!
!  Output:
!    File      : Next line in the file contains the data connected with the keyword KeyWd.
!
!  Calls:
!    MulaRdNlst
!    Normalize

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nUnit
character(len=*), intent(in) :: KeyWd
logical(kind=iwp), intent(in) :: rwnd
logical(kind=iwp), intent(out) :: Exists
integer(kind=iwp) :: BlInd, i, istatus, KLen
character(len=80) :: InLine, OutLine
character(len=32) :: TmpKey

InLine = KeyWd
call Normalize(InLine,OutLine)
! Index of first blank character from the end:
BlInd = 81
do i=80,1,-1
  if (Outline(i:i) /= ' ') exit
  BlInd = i
end do
if (BlInd > 80) then
  write(u6,*) 'KEYWORD: KeyWd is too long.'
  write(u6,*) 'KeyWd:',KeyWd
  write(u6,*) 'After normalization (OutLine):'
  write(u6,*) OutLine
  call abend()
end if
Exists = .false.
if (BlInd == 1) return

KLen = BlInd-1
TmpKey = OutLine(1:KLen)

if (rwnd) then
  rewind(nUnit)
  call MulaRdNlst(nUnit,'MULA')
end if

read(nUnit,'(a80)',iostat=istatus) InLine
if (istatus == 0) then
  call Normalize(InLine,OutLine)
  do while (OutLine(1:KLen) /= TmpKey(1:KLen))
    read(nUnit,'(a80)',iostat=istatus) InLine
    if (istatus /= 0) exit
    call Normalize(InLine,OutLine)
  end do
end if
if (istatus > 0) then
  write(u6,*) ' I/O error on unit nUnit=',nUnit
  call abend()
end if
Exists = OutLine(1:KLen) == TmpKey(1:KLen)

return

end subroutine KeyWord

function StrToDble(InString)
!  Purpose:
!    Convert a number in string InString to Real.
!
!  Input:
!    InString   : Character string.
!
!  Output:
!    xNum       : Real.

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: StrToDble
character(len=*), intent(in) :: InString
integer(kind=iwp) :: length
real(kind=wp) :: xNum
character(len=2) :: TwoDigits
character(len=1) :: OneDigit

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

function iStrToInt(InLine)
!  Purpose:
!    Convert a number in string InLine to integer.
!
!  Input:
!    InLine   : Character string.
!
!  Output:
!    num      : Integer.

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iStrToInt
character(len=*), intent(in) :: InLine
integer(kind=iwp) :: i, isum, length, nPos, num
logical(kind=iwp) :: minus
character :: ch

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
    isum = isum+(index('0123456789',ch)-1)*10**(i)
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

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: line
character(len=*), intent(out) :: line2
integer(kind=iwp), parameter :: ntrans = 256
integer(kind=iwp) :: i, icalled = 0, ista, itrans(ntrans), j, len1, len2, len2max
logical(kind=iwp) :: wrdend
character :: ch
character(len=*), parameter :: lcase = 'abcdefghijklmnopqrstuvwxyz', &
                               ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

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
    exit
  end if
end do
line2 = ' '
if (ista == 0) return
wrdend = .false.
len2 = 0
do i=ista,len1
  ch = line(i:i)
  if (ch /= ' ') then
    if (wrdend) then
      if (len2 > len2max-1) then
        write(u6,*) ' WARNING: Line truncated in NORMALIZE.'
        exit
      end if
      line2(len2+1:len2+1) = ' '
      len2 = len2+1
    end if
    if (len2 > len2max-1) then
      write(u6,*) ' WARNING: Line truncated in NORMALIZE.'
      exit
    end if
    ch = char(itrans(ichar(ch)))
    line2(len2+1:len2+1) = ch
    len2 = len2+1
    wrdend = .false.
  else
    wrdend = .true.
  end if
end do

return

end subroutine Normalize

!end module IOTools
