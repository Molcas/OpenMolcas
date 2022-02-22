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
! Copyright (C) 1991, Bjorn O. Roos                                    *
!               1991, Roland Lindh                                     *
!               Valera Veryazov                                        *
!***********************************************************************

function Lbl2Nr(Atom)

use Isotopes, only: MaxAtomNum, PTab
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Lbl2Nr
character(len=*), intent(in) :: Atom
integer(kind=iwp) :: i, lAtom
character(len=2) :: TmpLbl(2)

if ((len_trim(Atom) > 2) .or. (len_trim(Atom) <= 0)) then
  call WarningMessage(2,'The atom label;-->'//Atom(1:4)//'<--; is not a proper string to define an element.')
  call Quit_OnUserError()
end if

Lbl2Nr = -1

lAtom = 0
do i=1,80
  if (Atom(i:i) == ' ') exit
  lAtom = lAtom+1
end do
if (lAtom > 2) lAtom = 2
TmpLbl(1) = '  '
if (lAtom > 0) then
  if (lAtom == 1) then
    TmpLbl(1)(2:2) = Atom(1:1)
  else
    TmpLbl(1) = Atom(1:2)
  end if
  call UpCase(TmpLbl(1))
  do i=0,MaxAtomNum
    TmpLbl(2) = PTab(i)
    call UpCase(TmpLbl(2))
    if (TmpLbl(2) == TmpLbl(1)) then
      Lbl2Nr = i
      exit
    end if
  end do
end if
if (Lbl2Nr == -1) then
  call WarningMessage(2,'The atom label;-->'//Atom(1:4)//'<--; does not define an element.')
  call Quit_OnUserError()
end if

return

end function Lbl2Nr
