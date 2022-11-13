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

subroutine Pr2D(xyz2d,nT,nRys,la,lb,lc,ld,IfGrad,iPrint)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, iPrint
real(kind=wp), intent(in) :: xyz2d(nT,nRys,0:la+1,0:lb+1,0:lc+1,0:ld+1,3)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp) :: ia, ib, ic, iCar, id, ja, jb, jc, jd
character(len=80) :: Label
character(len=*), parameter :: ch(3) = [',x)',',y)',',z)']
real(kind=wp), external :: DDot_

write(u6,*)
write(u6,*) ' Printing the 2d-integrals'
write(u6,*)

Label = ' '
ja = 0
if (IfGrad(1,1) .or. IfGrad(2,1) .or. IfGrad(3,1)) ja = 1
jb = 0
if (IfGrad(1,2) .or. IfGrad(2,2) .or. IfGrad(3,1)) jb = 1
jc = 0
if (IfGrad(1,3) .or. IfGrad(2,3) .or. IfGrad(3,3)) jc = 1
jd = 0
if (IfGrad(1,4) .or. IfGrad(2,4) .or. IfGrad(3,4)) jd = 1
do ia=0,la+ja
  if (ia > la) jb = 0
  do ib=0,lb+jb
    if ((ia > la) .or. (ib > lb)) jc = 0
    do ic=0,lc+jc
      do id=0,ld+jd
        do iCar=1,3
          if ((ja == 1) .and. (ia == la+ja) .and. (.not. IfGrad(iCar,1))) cycle
          if ((jb == 1) .and. (ib == lb+jb) .and. (.not. IfGrad(iCar,2))) cycle
          if ((jc == 1) .and. (ic == lc+jc) .and. (.not. IfGrad(iCar,3))) cycle
          if ((jd == 1) .and. (id == ld+jd) .and. (.not. IfGrad(iCar,4))) cycle
          write(Label,'(A,4(I1,A))') ' xyz2D0(',ia,',',ib,',',ic,',',id,ch(iCar)
          if (iPrint >= 99) then
            call RecPrt(Label,' ',xyz2d(:,:,ia,ib,ic,id,iCar),nT,nRys)
          else
            write(u6,'(A)') Label
            write(u6,*) DDot_(nT*nRys,xyz2d(:,:,ia,ib,ic,id,iCar),1,xyz2d(:,:,ia,ib,ic,id,iCar),1)
          end if
        end do
      end do
    end do
  end do
end do

return

end subroutine Pr2D
