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

subroutine Lowdin_LP(S,C,nDim)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(in) :: S(nDim,nDim)
real(kind=wp), intent(out) :: C(nDim,nDim)
integer(kind=iwp) :: i, ioff, ip_B, ip_S, j
#include "WrkSpc.fh"

call Allocate_Work(ip_S,nDim*(nDim+1)/2)
call Allocate_Work(ip_B,nDim**2)

do i=1,nDim
  do j=1,i
    ioff = ip_S-1+i*(i-1)/2+j
    Work(iOff) = S(i,j)
  end do
end do
call dcopy_(nDim**2,[Zero],0,Work(ip_B),1)
call dcopy_(nDim,[One],0,Work(ip_B),nDim+1)

call Lowdin_LP_(Work(ip_S),C,nDim,Work(ip_B))

call Free_Work(ip_B)
call Free_Work(ip_S)

return

end subroutine Lowdin_LP
