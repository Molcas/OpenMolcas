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

subroutine Pr2D(xyz2D,nT,nRys,la,lb,lc,ld,IfGrad,iPrint)

implicit real*8(a-h,o-z)
real*8 xyz2d(nT,nRys,0:la+1,0:lb+1,0:lc+1,0:ld+1,3)
logical IfGrad(3,4)
character Label*80, ch(3)*3
data ch/',x)',',y)',',z)'/

write(6,*)
write(6,*) ' Printing the 2d-integrals'
write(6,*)

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
          if ((ja == 1) .and. (ia == la+ja) .and. (.not. IfGrad(iCar,1))) Go To 51
          if ((jb == 1) .and. (ib == lb+jb) .and. (.not. IfGrad(iCar,2))) Go To 51
          if ((jc == 1) .and. (ic == lc+jc) .and. (.not. IfGrad(iCar,3))) Go To 51
          if ((jd == 1) .and. (id == ld+jd) .and. (.not. IfGrad(iCar,4))) Go To 51
          write(Label,'(A,4(I1,A))') ' xyz2D0(',ia,',',ib,',',ic,',',id,ch(iCar)
          if (iPrint >= 99) then
            call RecPrt(Label,' ',xyz2d(1,1,ia,ib,ic,id,iCar),nT,nRys)
          else
            write(6,'(A)') Label
            write(6,*) DDot_(nT*nRys,xyz2d(1,1,ia,ib,ic,id,iCar),1,xyz2d(1,1,ia,ib,ic,id,iCar),1)
          end if
51        continue
        end do
      end do
    end do
  end do
end do

return

end subroutine Pr2D
