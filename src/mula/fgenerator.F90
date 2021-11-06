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

subroutine FGenerator(nMat,F,iCre,iAnn,trgrd,m_ord,mx_ord,nosc)

integer nMat(0:m_Ord,nOsc)
real*8 trgrd(3,nOsc)
real*8 F(0:m_ord,0:mx_ord,3)
real*8 sqr(0:50)
integer iCre(0:m_ord,nosc)
integer iAnn(0:m_ord,nosc)

!F = 0.0d0
call dcopy_((m_ord+1)*(mx_ord+1)*3,[0.0d0],0,F,1)
do i=0,50
  sqr(i) = sqrt(dble(i)/2.0d0)
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
