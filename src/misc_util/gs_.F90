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

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nInter, nVec
real(kind=wp) :: T(nInter,nVec), Thr
integer(kind=iwp) :: i, j
real(kind=wp) :: XX, XY
real(kind=wp), external :: DDot_

do i=1,nVec

  ! Order the vectors according to the diagonal value.

  call GS_Order(T(1,i),nInter,nVec-i+1)

  ! Normalize the vector

  XX = sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
  !write(u6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    call DScal_(nInter,One/XX,T(1,i),1)
  else
    call FZero(T(1,i),nInter)
#   ifdef _DEBUGPRINT_
    call RecPrt('GS_: T',' ',T,nInter,nInter)
#   endif
    cycle
  end if

  ! Orthogonalize against the previous vectors

  ! |Y(new)>=|Y> - <X|Y> * |X>

  do j=1,i-1
    XY = DDot_(nInter,T(1,i),1,T(1,j),1)
    !if (abs(XY) > Thr) then
    !  write(u6,*) 'GS_: j,XY=',j,XY
    call DaXpY_(nInter,-XY,T(1,j),1,T(1,i),1)
    !  XY = DDot_(nInter,T(1,i),1,T(1,j),1)
    !  write(u6,*) 'GS_: j,XY=',j,XY
    !end if
  end do

  ! Renormalize

  XX = sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
  !write(u6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    call DScal_(nInter,One/XX,T(1,i),1)
  else
    call FZero(T(1,i),nInter)
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('GS_: T',' ',T,nInter,nInter)
# endif
end do

return

end subroutine GS_
