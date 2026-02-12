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

subroutine RdSupS(LuInput,n,iBuff)

use output_ras, only: LF

implicit none
!***********************************************************************
!                                                                      *
!     Purpose:                                                         *
!     Read supersymmetry input.                                        *
!                                                                      *
!***********************************************************************
integer LuInput, n
integer iBuff(*)
integer is(288), ie(288)
character(len=288) Line
integer K, I, nRepeat, iLast, m, l, iZ

!----------------------------------------------------------------------*
!     Start procedure, initialize data counter                         *
!----------------------------------------------------------------------*
n = 0
k = -1
!---  Read next line as a chatacter string  ---------------------------*
100 read(LuInput,'(A)',end=900) Line
!---  Left adjust line  -----------------------------------------------*
Line = adjustl(Line)
if (Line(1:1) == ' ') goto 100
if (Line(1:1) == '*') goto 100
!---  Remove multiple intervening blanks  -----------------------------*
do i=1,287
  nRepeat = 0
  do while ((Line(i:i+1) == '  ') .and. (nRepeat < 288))
    nRepeat = nRepeat+1
    Line(i:287) = Line(i+1:288)
    Line(288:288) = ' '
  end do
end do
!---  Insert commas as the only valid separators  ---------------------*
do i=2,287
  if (Line(i:i) == ' ') then
    if ((Line(i-1:i-1) /= ' ') .and. (Line(i-1:i-1) /= ',')) then
      if ((Line(i+1:i+1) /= ' ') .and. (Line(i+1:i+1) /= ',')) Line(i:i) = ','
    end if
  end if
end do
!---  Get the last noblank character  ---------------------------------*
iLast = 0
do i=1,288
  if (Line(i:i) /= ' ') iLast = i
end do
!---  Initialize markers  ---------------------------------------------*
do i=1,288
  is(i) = 0
  ie(i) = 0
end do
!---  Divide the line into substrings  --------------------------------*
m = 1
is(m) = 0
do i=1,iLast
  if (Line(i:i) == ',') then
    m = m+1
    is(m) = i
  end if
end do
m = 0
do i=1,iLast
  if (Line(i:i) == ',') then
    m = m+1
    ie(m) = i
  end if
end do
if (Line(iLast:iLast) /= ',') m = m+1
ie(m) = iLast+1
!---  Read by substrings  ---------------------------------------------*
do i=1,m
  l = ie(i)-is(i)
  iz = 0
  if (l > 2) then
    read(Line(is(i)+1:ie(i)-1),*,err=910) iz
  else if ((l == 2) .and. (Line(is(i)+1:ie(i)-1) /= ' ')) then
    read(Line(is(i)+1:ie(i)-1),*,err=910) iz
  end if
  if (k == -1) then
    k = k+1
    n = iz
  else if (k < n) then
    k = k+1
    iBuff(k) = iz
  end if
end do
!---  If necessary continue by next line  -----------------------------*
if (k < n) goto 100
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
return

!----------------------------------------------------------------------*
!     Error exit                                                       *
!----------------------------------------------------------------------*
900 write(LF,*)
write(LF,'(6X,A)') ' RASSCF was reading supersymmetry input from'
write(LF,'(6X,A)') 'the input file, when some error occurred,'
write(LF,'(6X,A)') 'probably end of file.'
write(LF,*)
call Quit_OnUserError()
910 write(LF,*)
write(LF,'(6X,A)') ' RASSCF was reading supersymmetry input from'
write(LF,'(6X,A)') 'the input file, when some error occurred.'
write(LF,'(6X,A)') 'Some of the input data seems to be in error.'
write(LF,*)
call Quit_OnUserError()

end subroutine RdSupS
