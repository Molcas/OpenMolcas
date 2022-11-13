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

subroutine GS_Order(T,nInter,nVec)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter, nVec
real(kind=wp), intent(inout) :: T(nInter,nVec)
integer(kind=iwp) :: i, iDiag, j
real(kind=wp) :: Diag, DiagMax
real(kind=wp), external :: DDot_

#ifdef _DEBUGPRINT_
call RecPrt('GS_Order: T(orig)','(12F6.2)',T,nInter,nVec)
#endif
do j=1,nVec-1
  DiagMax = DDot_(nInter,T(:,j),1,T(:,j),1)
  iDiag = j
  do i=j+1,nVec
    Diag = DDot_(nInter,T(:,i),1,T(:,i),1)
    if (T(i,i) > DiagMax) then
      DiagMax = Diag
      iDiag = i
    end if
  end do
  if (iDiag /= j) call dswap_(nInter,T(:,iDiag),1,T(:,j),1)
end do
#ifdef _DEBUGPRINT_
call RecPrt('GS_Order: T(ordered)','(12F6.2)',T,nInter,nVec)
#endif

return

end subroutine GS_Order
