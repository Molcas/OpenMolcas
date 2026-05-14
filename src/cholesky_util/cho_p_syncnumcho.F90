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

subroutine Cho_P_SyncNumCho(NumCho,nSym)
!
! Purpose: sync global NumCho_G vector counter. On entry, NumCho is
!          the local counter (unchanged).

use Cholesky, only: Cho_Real_Par, NumCho_G, NumChT_G, tMisc
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, NumCho(nSym)
real(kind=wp) :: c1, c2, w1, w2

if (Cho_Real_Par) then
  call CWTime(c1,w1)
  NumCho_G(1:nSym) = NumCho(1:nSym)
  call Cho_GAIGop(NumCho_G,nSym,'max')
  NumChT_G = sum(NumCho_G(1:nSym))
  call CWTime(c2,w2)
  tMisc(1,5) = tMisc(1,5)+c2-c1
  tMisc(2,5) = tMisc(2,5)+w2-w1
end if

end subroutine Cho_P_SyncNumCho
