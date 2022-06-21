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

subroutine Util1(Alpha,Beta,nZeta,final,la,lb,Slaplb,Slamlb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the electric field integrals from                *
!         derivative integrals of the electric potential.              *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,3,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), Slaplb(nZeta,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2), &
       Slamlb(nZeta,(la)*(la+1)/2,(lb+1)*(lb+2)/2), Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2), &
       Slalbm(nZeta,(la+1)*(la+2)/2,(lb)*(lb+1)/2), Alpha(nZeta), Beta(nZeta)
character*80 Label
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*) ' In Util1 la,lb=',la,lb
  call RecPrt('Alpha',' ',Alpha,nZeta,1)
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ib=1,nElem(lb)
    write(Label,'(A,I2,A)') ' Slaplb(la,',ib,')'
    call RecPrt(Label,' ',Slaplb(1,1,ib),nZeta,nElem(la+1))
  end do
  if (la > 0) then
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A)') ' Slamlb(la,',ib,')'
      call RecPrt(Label,' ',Slamlb(1,1,ib),nZeta,nElem(la-1))
    end do
  end if
  do ib=1,nElem(lb+1)
    write(Label,'(A,I2,A)') ' Slalbp(la,',ib,')'
    call RecPrt(Label,' ',Slalbp(1,1,ib),nZeta,nElem(la))
  end do
  if (lb > 0) then
    do ib=1,nElem(lb-1)
      write(Label,'(A,I2,A)') ' Slalbm(la,',ib,')'
      call RecPrt(Label,' ',Slalbm(1,1,ib),nZeta,nElem(la))
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

        if ((ixa == 0) .and. (ixb == 0)) then

          do iZeta=1,nZeta
            final(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
          end do

        else if (ixa == 0) then

          do iZeta=1,nZeta
            final(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))- &
                                     dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
          end do

        else if (ixb == 0) then

          do iZeta=1,nZeta
            final(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)- &
                                     dble(ixa)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
          end do

        else

          do iZeta=1,nZeta
            final(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)- &
                                     dble(ixa)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))- &
                                     dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
          end do

        end if

        if ((iya == 0) .and. (iyb == 0)) then

          do iZeta=1,nZeta
            final(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
          end do

        else if (iya == 0) then

          do iZeta=1,nZeta
            final(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))- &
                                     dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
          end do

        else if (iyb == 0) then

          do iZeta=1,nZeta
            final(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)- &
                                     dble(iya)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
          end do

        else

          do iZeta=1,nZeta
            final(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)- &
                                     dble(iya)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))- &
                                     dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
          end do

        end if

        if ((iza == 0) .and. (izb == 0)) then

          do iZeta=1,nZeta
            final(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
          end do

        else if (iza == 0) then

          do iZeta=1,nZeta
            final(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))- &
                                     dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
          end do

        else if (izb == 0) then

          do iZeta=1,nZeta
            final(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)- &
                                     dble(iza)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
          end do

        else

          do iZeta=1,nZeta
            final(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)- &
                                     dble(iza)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)+ &
                                     Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))- &
                                     dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
          end do

        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(6,*) ' In Util1 la,lb=',la,lb
  do iElem=1,nElem(la)
    do jElem=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,') '
      call RecPrt(Label,' ',final(1,1,iElem,jElem),nZeta,3)
    end do
  end do
end if

return

end subroutine Util1
