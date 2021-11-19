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

subroutine FGenerator(nMat,F,iCre,iAnn,trgrd,m_Ord,mx_Ord,nOsc)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m_Ord, nOsc, nMat(0:m_Ord,nOsc), iCre(0:m_Ord,nOsc), iAnn(0:m_Ord,nOsc), mx_Ord
real(kind=wp), intent(out) :: F(0:m_Ord,0:mx_Ord,3)
real(kind=wp), intent(in) :: trgrd(3,nOsc)
integer(kind=iwp) :: i, iCar, iOrd, iOsc, jOrd
real(kind=wp) :: sqr(0:50)

F(:,:,:) = Zero
do i=0,50
  sqr(i) = sqrt(real(i,kind=wp)*Half)
end do
do iCar=1,3
  do iOsc=1,nOsc
    do iOrd=1,m_Ord
      jOrd = iAnn(iOrd,iOsc)
      if (jOrd >= 0) then
        F(iOrd,jOrd,iCar) = F(iOrd,jOrd,iCar)+sqr(nMat(iOrd,iOsc))*trgrd(iCar,iOsc)
      end if
    end do
  end do
  do iOsc=1,nOsc
    do iOrd=0,m_Ord
      jOrd = iCre(iord,iosc)
      if (jOrd >= 0) then
        F(iOrd,jOrd,iCar) = F(iOrd,jOrd,iCar)+sqr(nMat(jOrd,iOsc))*trgrd(iCar,iOsc)
      end if
    end do
  end do
end do

end subroutine FGenerator
