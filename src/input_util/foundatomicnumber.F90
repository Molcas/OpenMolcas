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

subroutine FoundAtomicNumber(LuWr,Symb,NAT,iErr)

use isotopes, only: MaxAtomNum, PTab
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LuWr
character(len=6), intent(inout) :: Symb
integer(kind=iwp), intent(out) :: NAT, iErr
integer(kind=iwp) :: i

if ((Symb(1:1) <= 'z') .and. (Symb(1:1) >= 'a')) Symb(1:1) = char(ichar(Symb(1:1))-32)
if ((Symb(2:2) <= 'Z') .and. (Symb(2:2) >= 'A')) Symb(2:2) = char(ichar(Symb(2:2))+32)
iErr = 1

do i=1,MaxAtomNum
  if (Symb(1:2) == adjustl(PTab(i))) then
    iErr = 0
    NAT = i
    return
  end if
end do
do i=1,MaxAtomNum
  if (' '//Symb(1:1) == PTab(i)) then
    iErr = 0
    NAT = i
    return
  end if
end do

! --- Z atoms are ghost atoms that will disappear in seward
if (Symb(1:1) == 'Z') then
  iErr = 0
  NAT = -1
  return
end if
! --- X atoms are dummy atoms that will appear in seward
if (Symb(1:1) == 'X') then
  iErr = 0
  NAT = 0
  return
end if

write(LuWr,*) '   [FoundAtomicNumber]: Wrong atomic symbol !'

return

end subroutine FoundAtomicNumber
