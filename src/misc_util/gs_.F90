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

subroutine GS_(T,nInter,nVec,Thr)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter, nVec
real(kind=wp), intent(inout) :: T(nInter,nVec)
real(kind=wp), intent(in) :: Thr
integer(kind=iwp) :: i, j
real(kind=wp) :: XX, XY
real(kind=wp), external :: DDot_

do i=1,nVec

  ! Order the vectors according to the diagonal value.

  call GS_Order(T(:,i:),nInter,nVec-i+1)

  ! Normalize the vector

  XX = sqrt(DDot_(nInter,T(:,i),1,T(:,i),1))
  !write(u6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    T(:,i) = T(:,i)/XX
  else
    T(:,i) = Zero
#   ifdef _DEBUGPRINT_
    call RecPrt('GS_: T',' ',T,nInter,nInter)
#   endif
    cycle
  end if

  ! Orthogonalize against the previous vectors

  ! |Y(new)>=|Y> - <X|Y> * |X>

  do j=1,i-1
    XY = DDot_(nInter,T(:,i),1,T(:,j),1)
    !if (abs(XY) > Thr) then
    !  write(u6,*) 'GS_: j,XY=',j,XY
    T(:,i) = T(:,i)-XY*T(:,j)
    !  XY = DDot_(nInter,T(1,i),1,T(1,j),1)
    !  write(u6,*) 'GS_: j,XY=',j,XY
    !end if
  end do

  ! Renormalize

  XX = sqrt(DDot_(nInter,T(:,i),1,T(:,i),1))
  !write(u6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    T(:,i) = T(:,i)/XX
  else
    T(:,i) = Zero
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('GS_: T',' ',T,nInter,nInter)
# endif
end do

return

end subroutine GS_
