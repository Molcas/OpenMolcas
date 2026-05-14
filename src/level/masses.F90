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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine MASSES(IAN,IMN,ANAME,MASS)
!***********************************************************************
!** For isotope with (input) atomic number IAN and mass number IMN,
!  return (output):  (i) as the right-adjusted 2-character variable ANAME
!  the alphabetic symbol for that element, and (ii) the atomic mass MASS
!  [amu].
!** If the input value of IMN does not equal one of the tabulated values
!  for atomic species IAN, return the abundance-averaged standard atomic
!  weight of that atom
!***********************************************************************

use Isotopes, only: ElementList, Initialize_Isotopes, MaxAtomNum
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ian
integer(kind=iwp), intent(inout) :: imn
character(len=2), intent(out) :: ANAME
real(kind=wp), intent(out) :: mass
integer(kind=iwp) :: i, n

call Initialize_Isotopes()

if ((IAN <= 0) .or. (IAN > MaxAtomNum)) then
  MASS = Zero
  ANAME = 'XX'
  IMN = 0
  write(u6,601) IAN
  return
else
  ANAME = adjustr(ElementList(IAN)%Symbol)
end if
if ((IAN == 1)) then
  ! Special case: insert common name for deuterium or tritium
  if (IMN == 2) ANAME = ' D'
  if (IMN == 3) ANAME = ' T'
end if
MASS = -One
do I=1,size(ElementList(IAN)%Isotopes)
  if (IMN == ElementList(IAN)%Isotopes(I)%A) then
    MASS = ElementList(IAN)%Isotopes(I)%m
    exit
  end if
end do
if (MASS < Zero) then
  n = ElementList(IAN)%Natural
  if (n < 1) then
    MASS = ElementList(IAN)%Isotopes(1)%m
  else
    MASS = sum(ElementList(IAN)%Isotopes(1:n)%m*ElementList(IAN)%Isotopes(1:n)%x)
  end if
  if (IMN /= 0) then
    write(u6,602) ANAME,IMN
    IMN = 0
  end if
end if

return

601 format(' *** Isotopes database does not include Atomic Number=',i4)
602 format(' *** Isotopes database does not include ',A2,'(',i3,'), so use average atomic mass.')

end subroutine MASSES
