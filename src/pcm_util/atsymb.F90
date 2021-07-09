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

function AtSymb(I)
! ATSYMB(I) is the Atomic Symbol corresponding to the Atomic number I
! if I=0  Atsymb=Bq
! if I=-1 Atsymb=X

use Definitions, only: iwp

use isotopes, only: PTab

implicit none
character(len=2) :: AtSymb
integer(kind=iwp), intent(in) :: I

if (I > 0) then
  AtSymb = PTab(i)
  return
end if
if (I == -1) AtSymb = ' X'
if (I == 0) AtSymb = 'Bq'

return

end function AtSymb
