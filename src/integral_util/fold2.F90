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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine Fold2(nSym,nBas,A,B)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(*)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp) :: iBas, iOff1, iOff2, iSym, mBas

iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  mBas = nBas(iSym)
  do iBas=1,mBas
    B(iOff2+1:iOff2+iBas) = A(iOff1+1:iOff1+iBas)
    iOff1 = iOff1+mBas
    iOff2 = iOff2+iBas
  end do
end do

return

end subroutine Fold2
