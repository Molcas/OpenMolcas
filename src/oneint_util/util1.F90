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

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Alpha(nZeta), Beta(nZeta), Slaplb(nZeta,nTri_Elem1(la+1),nTri_Elem1(lb)), &
                             Slamlb(nZeta,nTri_Elem1(la-1),nTri_Elem1(lb)), Slalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1)), &
                             Slalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1))
real(kind=wp), intent(out) :: rFinal(nZeta,3,nTri_Elem1(la),nTri_Elem1(lb))
#include "print.fh"
integer(kind=iwp) :: ib, iElem, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, jElem
character(len=80) :: Label

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In Util1 la,lb=',la,lb
  call RecPrt('Alpha',' ',Alpha,nZeta,1)
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ib=1,nTri_Elem1(lb)
    write(Label,'(A,I2,A)') ' Slaplb(la,',ib,')'
    call RecPrt(Label,' ',Slaplb(:,:,ib),nZeta,nTri_Elem1(la+1))
  end do
  if (la > 0) then
    do ib=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A)') ' Slamlb(la,',ib,')'
      call RecPrt(Label,' ',Slamlb(:,:,ib),nZeta,nTri_Elem1(la-1))
    end do
  end if
  do ib=1,nTri_Elem1(lb+1)
    write(Label,'(A,I2,A)') ' Slalbp(la,',ib,')'
    call RecPrt(Label,' ',Slalbp(:,:,ib),nZeta,nTri_Elem1(la))
  end do
  if (lb > 0) then
    do ib=1,nTri_Elem1(lb-1)
      write(Label,'(A,I2,A)') ' Slalbm(la,',ib,')'
      call RecPrt(Label,' ',Slalbm(:,:,ib),nZeta,nTri_Elem1(la))
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

        if ((ixa /= 0) .and. (ixb /= 0)) then
          rFinal(:,1,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa+1,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))- &
                                real(ixa,kind=wp)*Slamlb(:,C_Ind(la-1,ixa-1,iza),ipb)- &
                                real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb))
        else if (ixa /= 0) then
          rFinal(:,1,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa+1,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))- &
                                real(ixa,kind=wp)*Slamlb(:,C_Ind(la-1,ixa-1,iza),ipb)
        else if (ixb /= 0) then
          rFinal(:,1,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa+1,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))- &
                                real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb))
        else
          rFinal(:,1,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa+1,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))
        end if

        if ((iya /= 0) .and. (iyb /= 0)) then
          rFinal(:,2,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))- &
                                real(iya,kind=wp)*Slamlb(:,C_Ind(la-1,ixa,iza),ipb)- &
                                real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb))
        else if (iya /= 0) then
          rFinal(:,2,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))- &
                                real(iya,kind=wp)*Slamlb(:,C_Ind(la-1,ixa,iza),ipb)
        else if (iyb /= 0) then
          rFinal(:,2,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))- &
                                real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb))
        else
          rFinal(:,2,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))
        end if

        if ((iza /= 0) .and. (izb /= 0)) then
          rFinal(:,3,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza+1),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))- &
                                real(iza,kind=wp)*Slamlb(:,C_Ind(la-1,ixa,iza-1),ipb)- &
                                real(izb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1))
        else if (iza /= 0) then
          rFinal(:,3,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza+1),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))- &
                                real(iza,kind=wp)*Slamlb(:,C_Ind(la-1,ixa,iza-1),ipb)
        else if (izb /= 0) then
          rFinal(:,3,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza+1),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))- &
                                real(izb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1))
        else
          rFinal(:,3,ipa,ipb) = Two*Alpha*Slaplb(:,C_Ind(la+1,ixa,iza+1),ipb)+Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))
        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(u6,*) ' In Util1 la,lb=',la,lb
  do iElem=1,nTri_Elem1(la)
    do jElem=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,') '
      call RecPrt(Label,' ',rFinal(:,:,iElem,jElem),nZeta,3)
    end do
  end do
end if

return

end subroutine Util1
