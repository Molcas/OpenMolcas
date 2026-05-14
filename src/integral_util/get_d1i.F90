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

subroutine Get_D1I(CMO,D1It,D1I,nish,nbas,nsym)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: D1It(*), D1I(*)
integer(kind=iwp), intent(in) :: nSym, nish(nsym), nbas(nsym)
integer(kind=iwp) :: i, iBas, iOff1, iOff2, iOrb, iSym, j, k
real(kind=wp) :: rSum

iOff1 = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  iOrb = nIsh(iSym)
  if (iBas /= 0) then
    iOff2 = iOff1
    do i=1,iBas
      do j=1,iBas
        rSum = Zero
        do k=0,iOrb-1
          rSum = rSum+Two*CMO(iOff1+k*iBas+i)*CMO(iOff1+k*iBas+j)
        end do
        D1I(iOff2+j) = rSum
      end do
      iOff2 = iOff2+iBas
    end do
    iOff1 = iOff1+iBas*iBas
  end if
end do
call Fold2(nsym,nBas,D1I,D1It)

return

end subroutine Get_D1I
