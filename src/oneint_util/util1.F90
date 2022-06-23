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

subroutine Util1(Alpha,Beta,nZeta,rFinal,la,lb,Slaplb,Slamlb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the electric field integrals from                *
!         derivative integrals of the electric potential.              *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nZeta, la, lb
real(kind=wp) :: Alpha(nZeta), Beta(nZeta), rFinal(nZeta,3,nTri_Elem1(la),nTri_Elem1(lb)), &
                 Slaplb(nZeta,nTri_Elem1(la+1),nTri_Elem1(lb)), Slamlb(nZeta,nTri_Elem1(la-1),nTri_Elem1(lb)), &
                 Slalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1)), Slalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1))
#include "print.fh"
integer(kind=iwp) :: ib, iElem, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, iZeta, jElem
character(len=80) :: Label
! Statement function for cartesian index
integer(kind=iwp) :: Ind, nElem, ixyz, ix, iz
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In Util1 la,lb=',la,lb
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
            rFinal(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
          end do

        else if (ixa == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))- &
                                      real(ixb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
          end do

        else if (ixb == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)- &
                                      real(ixa,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
          end do

        else

          do iZeta=1,nZeta
            rFinal(iZeta,1,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa+1,iza),ipb)- &
                                      real(ixa,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa-1,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))- &
                                      real(ixb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
          end do

        end if

        if ((iya == 0) .and. (iyb == 0)) then

          do iZeta=1,nZeta
            rFinal(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
          end do

        else if (iya == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))- &
                                      real(iyb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
          end do

        else if (iyb == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)- &
                                      real(iya,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
          end do

        else

          do iZeta=1,nZeta
            rFinal(iZeta,2,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza),ipb)- &
                                      real(iya,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa,iza),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))- &
                                      real(iyb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
          end do

        end if

        if ((iza == 0) .and. (izb == 0)) then

          do iZeta=1,nZeta
            rFinal(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
          end do

        else if (iza == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))- &
                                      real(izb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
          end do

        else if (izb == 0) then

          do iZeta=1,nZeta
            rFinal(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)- &
                                      real(iza,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
          end do

        else

          do iZeta=1,nZeta
            rFinal(iZeta,3,ipa,ipb) = Two*Alpha(iZeta)*Slaplb(iZeta,Ind(la+1,ixa,iza+1),ipb)- &
                                      real(iza,kind=wp)*Slamlb(iZeta,Ind(la-1,ixa,iza-1),ipb)+ &
                                      Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))- &
                                      real(izb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
          end do

        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(u6,*) ' In Util1 la,lb=',la,lb
  do iElem=1,nElem(la)
    do jElem=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,') '
      call RecPrt(Label,' ',rFinal(1,1,iElem,jElem),nZeta,3)
    end do
  end do
end if

return

end subroutine Util1
