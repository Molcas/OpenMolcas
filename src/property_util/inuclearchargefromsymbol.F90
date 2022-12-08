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
! Copyright (C) 2000, Per-Olof Widmark                                 *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

function iNuclearChargeFromSymbol(Symbol)
!***********************************************************************
!                                                                      *
! This function returns the nuclear charge of an atom based on the     *
! chemical symbol.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
! Written: Feb. 2000                                                   *
! Modified: March 2017, Ignacio Fdez. Galvan (use PTab)                *
!                                                                      *
!***********************************************************************

use Isotopes, only: MaxAtomNum, PTab
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iNuclearChargeFromSymbol
character(len=*), intent(in) :: Symbol
integer(kind=iwp) :: i, idx
character(len=2) :: Sym1, Sym2

!----------------------------------------------------------------------*
! Locate symbol in table.                                              *
!----------------------------------------------------------------------*
idx = 0
Sym1 = adjustl(Symbol)
call UpCase(Sym1)
do i=1,MaxAtomNum
  Sym2 = adjustl(PTab(i))
  call UpCase(Sym2)
  if (Sym1 == Sym2) idx = i
end do
!----------------------------------------------------------------------*
! Are we successful?                                                   *
!----------------------------------------------------------------------*
if (idx == 0) then
  write(u6,'(a)') '***'
  write(u6,'(a)') '*** iNuclearChargeFromSymbol: warning'
  write(u6,'(2a)') '***    unknown atom: ',Symbol
  write(u6,'(a)') '***'
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
iNuclearChargeFromSymbol = idx

return

end function iNuclearChargeFromSymbol
