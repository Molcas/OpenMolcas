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

implicit real*8(a-h,o-z)
real*8 T(nInter,nVec)
#include "real.fh"

do i=1,nVec

  ! Order the vectors according to the diagonal value.

  call GS_Order(T(1,i),nInter,nVec-i+1)

  ! Normalize the vector

  XX = sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
  !write(6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    call DScal_(nInter,One/XX,T(1,i),1)
  else
    call FZero(T(1,i),nInter)
    Go To 100
  end if

  ! Orthogonalize against the previous vectors

  ! |Y(new)>=|Y> - <X|Y> * |X>

  do j=1,i-1
    XY = DDot_(nInter,T(1,i),1,T(1,j),1)
    !if (abs(XY) > Thr) then
    !  write(6,*) 'GS_: j,XY=',j,XY
    call DaXpY_(nInter,-XY,T(1,j),1,T(1,i),1)
    !  XY = DDot_(nInter,T(1,i),1,T(1,j),1)
    !  write(6,*) 'GS_: j,XY=',j,XY
    !end if
  end do

  ! Renormalize

  XX = sqrt(DDot_(nInter,T(1,i),1,T(1,i),1))
  !write(6,*) 'GS_: i,XX=',i,XX
  if (XX > Thr) then
    call DScal_(nInter,One/XX,T(1,i),1)
  else
    call FZero(T(1,i),nInter)
  end if
100 continue
# ifdef _DEBUGPRINT_
  call RecPrt('GS_: T',' ',T,nInter,nInter)
# endif
end do

return

end subroutine GS_
