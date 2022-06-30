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
!               2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Util5(Beta,nZeta,rFinal,la,lb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the velocity quadrupole integrals from the       *
!         derivative integrals and dipole integrals.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!     Adapted from util2: Ignacio Fdez. Galvan, July 2019              *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Beta(nZeta), Slalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1),3), &
                             Slalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1),3)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iElem, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, jElem
character(len=80) :: Label

iRout = 211
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In Util5 la,lb=',la,lb
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'x)'
      call RecPrt(Label,' ',Slalbp(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'y)'
      call RecPrt(Label,' ',Slalbp(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'z)'
      call RecPrt(Label,' ',Slalbp(:,ia,ib,3),nZeta,1)
    end do
  end do
  if (lb > 0) then
    do ia=1,nTri_Elem1(la)
      do ib=1,nTri_Elem1(lb-1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'x)'
        call RecPrt(Label,' ',Slalbm(:,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'y)'
        call RecPrt(Label,' ',Slalbm(:,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'z)'
        call RecPrt(Label,' ',Slalbm(:,ia,ib,3),nZeta,1)
      end do
    end do
  end if
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = C_Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)

        rFinal(:,ipa,ipb,1) = Two*Beta*Two*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb),1)
        rFinal(:,ipa,ipb,2) = Two*Beta*(Slalbp(:,ipa,C_Ind(lb+1,ixb,izb),1)+Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb),2))
        rFinal(:,ipa,ipb,3) = Two*Beta*(Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1),1)+Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb),3))
        rFinal(:,ipa,ipb,4) = Two*Beta*Two*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb),2)
        rFinal(:,ipa,ipb,5) = Two*Beta*(Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1),2)+Slalbp(:,ipa,C_Ind(lb+1,ixb,izb),3))
        rFinal(:,ipa,ipb,6) = Two*Beta*Two*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1),3)

        if (ixb > 0) then
          rFinal(:,ipa,ipb,1) = rFinal(:,ipa,ipb,1)-real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb),1)*Two
          rFinal(:,ipa,ipb,2) = rFinal(:,ipa,ipb,2)-real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb),2)
          rFinal(:,ipa,ipb,3) = rFinal(:,ipa,ipb,3)-real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb),3)
        end if

        if (iyb > 0) then
          rFinal(:,ipa,ipb,2) = rFinal(:,ipa,ipb,2)-real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb),1)
          rFinal(:,ipa,ipb,4) = rFinal(:,ipa,ipb,4)-real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb),2)*Two
          rFinal(:,ipa,ipb,5) = rFinal(:,ipa,ipb,5)-real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb),3)
        end if

        if (izb > 0) then
          rFinal(:,ipa,ipb,3) = rFinal(:,ipa,ipb,3)-real(izb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1),1)
          rFinal(:,ipa,ipb,5) = rFinal(:,ipa,ipb,5)-real(izb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1),2)
          rFinal(:,ipa,ipb,6) = rFinal(:,ipa,ipb,6)-real(izb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1),3)*Two
        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(u6,*) ' In Util5 la,lb=',la,lb
  do iElem=1,nTri_Elem1(la)
    do jElem=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',xx) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',xy) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',xz) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',yy) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',yz) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,',zz) '
      call RecPrt(Label,' ',rFinal(:,iElem,jElem,6),nZeta,1)
    end do
  end do
end if

return

end subroutine Util5
