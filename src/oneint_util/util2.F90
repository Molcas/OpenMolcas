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

subroutine Util2(Beta,nZeta,final,la,lb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the orbital angular momentum integrals from the  *
!         derivative integrals and dipole integrals.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3), Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3), &
       Slalbm(nZeta,(la+1)*(la+2)/2,lb*(lb+1)/2,3), Beta(nZeta)
character*80 Label
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 211
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*) ' In Util2 la,lb=',la,lb
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ia=1,nElem(la)
    do ib=1,nElem(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'x)'
      call RecPrt(Label,' ',Slalbp(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'y)'
      call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'z)'
      call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
    end do
  end do
  if (lb > 0) then
    do ia=1,nElem(la)
      do ib=1,nElem(lb-1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'x)'
        call RecPrt(Label,' ',Slalbm(1,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'y)'
        call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'z)'
        call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
      end do
    end do
  end if
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
          final(iZeta,ipa,ipb,1) = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3))
          final(iZeta,ipa,ipb,2) = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1))
          final(iZeta,ipa,ipb,3) = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),1)-Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2))
        end do

        if (ixb > 0) then
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,2) = final(iZeta,ipa,ipb,2)-dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3)
            final(iZeta,ipa,ipb,3) = final(iZeta,ipa,ipb,3)+dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2)
          end do
        end if

        if (iyb > 0) then
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,3) = final(iZeta,ipa,ipb,3)-dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),1)
            final(iZeta,ipa,ipb,1) = final(iZeta,ipa,ipb,1)+dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),3)
          end do
        end if

        if (izb > 0) then
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,1) = final(iZeta,ipa,ipb,1)-dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
            final(iZeta,ipa,ipb,2) = final(iZeta,ipa,ipb,2)+dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1)
          end do
        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(6,*) ' In Util2 la,lb=',la,lb
  do iElem=1,nElem(la)
    do jElem=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',x) '
      call RecPrt(Label,' ',final(1,iElem,jElem,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',y) '
      call RecPrt(Label,' ',final(1,iElem,jElem,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',z) '
      call RecPrt(Label,' ',final(1,iElem,jElem,3),nZeta,1)
    end do
  end do
end if

return

end subroutine Util2
