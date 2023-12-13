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
! Copyright (C) 2014,2015, Ignacio Fdez. Galvan                        *
!***********************************************************************
!  Preprocess_UDC
!
!> @brief
!>   Quickly read the constraints to take some decisions
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Process the constraints in the \p Lu file, but just to find out whether
!> some constraints are present and their values.
!> This is used for:
!>   - Detecting if the "EDIFF" constraint is present and with a zero value,
!>     in which case a conical intersection algorithm could be activated
!>     (depending on the symmetry and spin in the runfile).
!>   - Detecting if any constraint is explicitly declared as "soft",
!>     this is recommended for some TS searches with numerical differentiation.
!>   - Detecting whether an explicit MEP/IRC constraint has been included,
!>     this will probably give problems.
!>
!> @param[in]  Lu     Unit number of the file with the constraints
!> @param[in]  iPrint Print level
!***********************************************************************

subroutine Preprocess_UDC(Lu,iPrint)

use Slapaf_Info, only: EDiffZero, iState, lSoft, MEP_Type, MEPCons, NADC, NADC
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu, iPrint
integer(kind=iwp) :: Error, i, iPos, j, nLines
real(kind=wp) :: EDiffValue
logical(kind=iwp) :: MECI_via_SLAPAF
character(len=180) :: EDiffName, Line1, Line2
character(len=180), external :: Get_Ln

EDiffName = ''
EDiffZero = .false.
iState(1) = 0
iState(2) = 0
MECI_via_SLAPAF = .false.

! An arbitrary initial value given to EDiffValue, this value is
! irrelevant since the only thing that matters is whether it is 0.0

EDiffValue = One

! Ugly hack to be able to restore the file to its original position
! (since Get_Ln may read more than one line)
! Count the lines until the end of the file and backspace

nLines = 0
do
  read(Lu,'(A)',IOSTAT=Error) Line1
  nLines = nLines+1
  if (Error /= 0) exit
end do
do i=1,nLines
  backspace(Lu)
end do

! Read the primitive constraints

! First read until the "VALUES" line
! Each line is split at the "=" sign by auto,
! so read two lines if there is no "=" sign

do
  Line1 = Get_Ln(Lu)
  call UpCase(Line1)
  Line1 = adjustl(Line1)
  if (Line1(1:4) == 'VALU') exit
  iPos = index(Line1,'=')
  if (iPos > 0) then
    Line2 = Line1(iPos+1:)
    Line1 = Line1(:iPos-1)
  else
    Line2 = Get_Ln(Lu)
    call UpCase(Line2)
  end if
  Line2 = adjustl(Line2)
  ! If a primitive is defined as "EDIFF", save its name
  if (Line2(1:4) == 'EDIF') then
    EDiffName = Line1
    iPos = index(Line2,' ')
    Line2 = Line2(iPos+1:)
    read(Line2,*,IOSTAT=Error) i
    if (Error /= 0) i = 0
    iPos = index(Line2,' ')
    Line2 = Line2(iPos+1:)
    read(Line2,*,IOSTAT=Error) j
    if (Error /= 0) j = 0
    iState(1) = max(i,j)
    iState(2) = min(i,j)
  end if
  ! If a primitive of the same type as that used by MEP/IRC
  ! is defined, signal it.
  if (Line2(1:4) == MEP_Type(1:4)) MEPCons = .true.
end do

! Now read the constraint values

! Read until the "END" line
! Each line is split at the "=" sign by auto,
! so read two lines if there is no "=" sign

do
  Line1 = Get_Ln(Lu)
  call UpCase(Line1)
  Line1 = adjustl(Line1)
  ! Read continuation lines
  do while (index(Line1,'&') /= 0)
    Line1 = Get_Ln(Lu)
    call UpCase(Line1)
    Line1 = adjustl(Line1)
  end do
  if (Line1(1:4) == 'END ') exit
  iPos = index(Line1,'=')
  if (iPos > 0) then
    Line2 = Line1(iPos+1:)
    Line1 = Line1(:iPos-1)
  else
    Line2 = Get_Ln(Lu)
    call UpCase(Line2)
  end if
  ! If the name matches the "EDIFF" primitive, read the value
  if (Line1 == EDiffName) then
    read(Line2,*,IOSTAT=Error) EDiffValue
    if (Error /= 0) EDiffValue = One
  end if
  ! If the constraint is explicitly "soft"
  if (index(Line2,'SOFT') /= 0) lSoft = .true.
end do

! Return the file to its original position
! (cannot use REWIND because this may be the full input for slapaf)
! First advance to the end of file and then backspace

do
  read(Lu,'(A)',IOSTAT=Error) Line1
  if (Error /= 0) exit
end do
do i=1,nLines
  backspace(Lu)
end do

! If an "EDIFF" constraint is being used with a value exactly 0.0,
! this may be a conical intersection search.
! If there is no "EDIFF" at all, there's no use computing NACs

if (EDiffValue == Zero) then
  EDiffZero = .true.
  MECI_via_SLAPAF = .true.
  call put_lscalar('MECI_via_SLAPAF ',MECI_via_SLAPAF)
  if (iPrint >= 6) then
    write(u6,*) 'Energy difference constraint with zero value.'
    write(u6,*) 'This may be a conical intersection search.'
  end if
else if (EDiffName(1:4) /= '    ') then
  if (iPrint >= 6) then
    write(u6,*) 'Energy difference constraint with non-zero value.'
    write(u6,*) 'This will not be a conical intersection search.'
  end if
else
  NADC = .false.
end if

end subroutine Preprocess_UDC
