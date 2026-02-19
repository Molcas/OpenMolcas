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
!***********************************************************************
!                                                                      *
!     Purpose:                                                         *
!     Read supersymmetry input.                                        *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: LuInput, n, iBuff(*)
integer(kind=iwp) :: I, ie(288), iLast, is(288), istatus, iZ, K, l, m, nRepeat
character(len=288) :: Line

!----------------------------------------------------------------------*
!     Start procedure, initialize data counter                         *
!----------------------------------------------------------------------*
n = 0
k = -1
!---  Read next line as a character string  ---------------------------*
do
  read(LuInput,'(A)',iostat=istatus) Line
  if (istatus < 0) call Error(istatus)
  !---  Left adjust line  ---------------------------------------------*
  Line = adjustl(Line)
  if ((Line(1:1) == ' ') .or.(Line(1:1) == '*')) cycle
  !---  Remove multiple intervening blanks  ---------------------------*
  do i=1,287
    nRepeat = 0
    do while ((Line(i:i+1) == '  ') .and. (nRepeat < 288))
      nRepeat = nRepeat+1
      Line(i:287) = Line(i+1:288)
      Line(288:288) = ' '
    end do
  end do
  !---  Insert commas as the only valid separators  -------------------*
  do i=2,287
    if (Line(i:i) == ' ') then
      if ((Line(i-1:i-1) /= ' ') .and. (Line(i-1:i-1) /= ',')) then
        if ((Line(i+1:i+1) /= ' ') .and. (Line(i+1:i+1) /= ',')) Line(i:i) = ','
      end if
    end if
  end do
  !---  Get the last noblank character  -------------------------------*
  iLast = 0
  do i=1,288
    if (Line(i:i) /= ' ') iLast = i
  end do
  !---  Initialize markers  -------------------------------------------*
  is(:) = 0
  ie(:) = 0
  !---  Divide the line into substrings  ------------------------------*
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
  !---  Read by substrings  -------------------------------------------*
  do i=1,m
    l = ie(i)-is(i)
    iz = 0
    if (l > 2) then
      read(Line(is(i)+1:ie(i)-1),*,iostat=istatus) iz
      if (istatus > 0) call Error(istatus)
    else if ((l == 2) .and. (Line(is(i)+1:ie(i)-1) /= ' ')) then
      read(Line(is(i)+1:ie(i)-1),*,iostat=istatus) iz
      if (istatus > 0) call Error(istatus)
    end if
    if (k == -1) then
      k = k+1
      n = iz
    else if (k < n) then
      k = k+1
      iBuff(k) = iz
    end if
  end do
  !---  If necessary continue by next line  ---------------------------*
  if (k >= n) exit
end do
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
return

contains

!----------------------------------------------------------------------*
!     Error exit                                                       *
!----------------------------------------------------------------------*
subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  write(u6,*)
  if (code < 0) then
    write(u6,'(6X,A)') ' RASSCF was reading supersymmetry input from'
    write(u6,'(6X,A)') 'the input file, when some error occurred,'
    write(u6,'(6X,A)') 'probably end of file.'
  else if (code > 0) then
    write(u6,'(6X,A)') ' RASSCF was reading supersymmetry input from'
    write(u6,'(6X,A)') 'the input file, when some error occurred.'
    write(u6,'(6X,A)') 'Some of the input data seems to be in error.'
  end if
  write(u6,*)
  call Quit_OnUserError()

end subroutine Error

end subroutine RdSupS
