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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Util4(nZeta,rFinal,la,lb,Elalbp,Elalb,Bcoor,Dcoor)
!***********************************************************************
!                                                                      *
! Object: to assemble the diamagnetic shielding integrals from         *
!         electric field integrals.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Half, c_in_au
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nZeta, la, lb
real(kind=wp) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),9), Elalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1),3), &
                 Elalb(nZeta,nTri_Elem1(la),nTri_Elem1(lb),3), Bcoor(3), Dcoor(3)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iComp, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, iZeta
real(kind=wp) :: BD(3), Fact, xCxD, xCyD, xCzD, yCxD, yCyD, yCzD, zCxD, zCyD, zCzD
character(len=80) :: Label
! Statement function for cartesian index
integer(kind=iwp) :: Ind, nElem, ixyz, ix, iz
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 231
iPrint = nPrint(iRout)

BD(1) = Bcoor(1)-Dcoor(1)
BD(2) = Bcoor(2)-Dcoor(2)
BD(3) = Bcoor(3)-Dcoor(3)
Fact = -1.0e-6_wp*Half/c_in_au**2
if (iPrint >= 99) then
  write(u6,*) ' In Util4 la,lb=',la,lb
  do ia=1,nElem(la)
    do ib=1,nElem(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,3),nZeta,1)
    end do
  end do
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Elalb(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalb(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalb(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalb(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalb(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalb(1,ia,ib,3),nZeta,1)
    end do
  end do
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        do iZeta=1,nZeta
          xCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1)+BD(1)*Elalb(iZeta,ipa,ipb,1)
          yCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2)+BD(1)*Elalb(iZeta,ipa,ipb,2)
          zCxD = Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3)+BD(1)*Elalb(iZeta,ipa,ipb,3)
          xCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),1)+BD(2)*Elalb(iZeta,ipa,ipb,1)
          yCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),2)+BD(2)*Elalb(iZeta,ipa,ipb,2)
          zCyD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3)+BD(2)*Elalb(iZeta,ipa,ipb,3)
          xCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1)+BD(3)*Elalb(iZeta,ipa,ipb,1)
          yCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)+BD(3)*Elalb(iZeta,ipa,ipb,2)
          zCzD = Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3)+BD(3)*Elalb(iZeta,ipa,ipb,3)
          rFinal(iZeta,ipa,ipb,1) = Fact*(+yCyD+zCzD)
          rFinal(iZeta,ipa,ipb,2) = Fact*(-yCxD)
          rFinal(iZeta,ipa,ipb,3) = Fact*(-zCxD)
          rFinal(iZeta,ipa,ipb,4) = Fact*(-xCyD)
          rFinal(iZeta,ipa,ipb,5) = Fact*(xCxD+zCzD)
          rFinal(iZeta,ipa,ipb,6) = Fact*(-zCyD)
          rFinal(iZeta,ipa,ipb,7) = Fact*(-xCzD)
          rFinal(iZeta,ipa,ipb,8) = Fact*(-yCzD)
          rFinal(iZeta,ipa,ipb,9) = Fact*(xCxD+yCyD)
        end do

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  do iComp=1,9
    write(Label,'(A,I2,A)') ' rFinal (',iComp,') '
    call RecPrt(Label,' ',rFinal(1,1,1,iComp),nZeta,nElem(la)*nELem(lb))
  end do
end if

return

end subroutine Util4
