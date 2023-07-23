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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nSym, NumCho(nSym)
#include "choglob.fh"
#include "cho_para_info.fh"
integer(kind=iwp) :: iSym
real(kind=wp) :: c1, c2, w1, w2

if (Cho_Real_Par) then
  call Cho_Timer(c1,w1)
  call iCopy(nSym,NumCho,1,NumCho_G,1)
  call Cho_GAIGop(NumCho_G,nSym,'max')
  NumChT_G = NumCho_G(1)
  do iSym=2,nSym
    NumChT_G = NumChT_G+NumCho_G(iSym)
  end do
  call Cho_Timer(c2,w2)
  call Cho_P_SyncNumCho_Time(c2-c1,w2-w1)
end if

end subroutine Cho_P_SyncNumCho
