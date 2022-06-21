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
! Copyright (C) 1991,2015, Roland Lindh                                *
!               2015, Lasse Kragh Soerensen                            *
!***********************************************************************

subroutine Util3(Beta,nZeta,final,la,lb,Slalbp,Slalb,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the orbital magnetic quadrupole integrals from   *
!         the derivative integrals and dipole integrals.               *
!                                                                      *
!     Author: Roland Lindh, Lasse Kragh Soerensen                      *
!             Dept. of Theoretical Chemistry,                          *
!             University of Uppsala, SWEDEN                            *
!             February '15                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,9), Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,6), &
       Slalb(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3), Slalbm(nZeta,(la+1)*(la+2)/2,lb*(lb+1)/2,6), Beta(nZeta)
! Notice CmbnMP has just 6 components instead of 9!!! (automatically assumes symmetry) Well fuck you CmbnMP
! This means Slalbp and Slalbm in reality only have 6 components....
! XX = 1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5 and ZZ=6
! There are only six components since zy d/di = yz d/di
! We still keep the 9 components in final
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character*80 Label
#endif
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
#ifdef _DEBUGPRINT_
nElem(ix) = (ix+1)*(ix+2)/2
#endif

#ifdef _DEBUGPRINT_
write(6,*) ' In Util3 la,lb=',la,lb
call RecPrt('Beta',' ',Beta,nZeta,1)
do ia=1,nElem(la)
  do ib=1,nElem(lb+1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'xx)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,1),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'xy)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'xz)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'yx)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'yy)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,4),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'yz)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,5),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'zx)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'zy)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,5),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalbp(',ia,',',ib,'zz)'
    call RecPrt(Label,' ',Slalbp(1,ia,ib,6),nZeta,1)
  end do
end do
do ia=1,nElem(la)
  do ib=1,nElem(lb)
    write(Label,'(A,I2,A,I2,A)') ' Slalb (',ia,',',ib,'x)'
    call RecPrt(Label,' ',Slalb(1,ia,ib,1),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalb (',ia,',',ib,'y)'
    call RecPrt(Label,' ',Slalb(1,ia,ib,2),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Slalb (',ia,',',ib,'z)'
    call RecPrt(Label,' ',Slalb(1,ia,ib,3),nZeta,1)
  end do
end do
if (lb > 0) then
  do ia=1,nElem(la)
    do ib=1,nElem(lb-1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'xx)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'xy)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'xz)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'yx)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'yy)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'yz)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'zx)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'zy)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Slalbm(',ia,',',ib,'zz)'
      call RecPrt(Label,' ',Slalbm(1,ia,ib,6),nZeta,1)
    end do
  end do
end if
#endif

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        do iZeta=1,nZeta
          temp_xx = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3))
          temp_xy = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),4)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),5))
          temp_xz = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),5)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),6))
          temp_yx = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1))
          temp_yy = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),5)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2))
          temp_yz = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),6)-Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3))
          temp_zx = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),1)-Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2))
          temp_zy = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),2)-Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),4))
          temp_zz = Two*Beta(iZeta)*(Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3)-Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),5))

          temp_x = Slalb(iZeta,ipa,Ind(lb,ixb,izb),1)
          temp_y = Slalb(iZeta,ipa,Ind(lb,ixb,izb),2)
          temp_z = Slalb(iZeta,ipa,Ind(lb,ixb,izb),3)

          ! xx term
          final(iZeta,ipa,ipb,1) = Two*temp_xx
          ! xy term
          final(iZeta,ipa,ipb,2) = Two*temp_xy-temp_z
          ! xz term
          final(iZeta,ipa,ipb,3) = Two*temp_xz+temp_y
          ! yx term
          final(iZeta,ipa,ipb,4) = Two*temp_yx+temp_z
          ! yy term
          final(iZeta,ipa,ipb,5) = Two*temp_yy
          ! yz term
          final(iZeta,ipa,ipb,6) = Two*temp_yz-temp_x
          ! zx term
          final(iZeta,ipa,ipb,7) = Two*temp_zx-temp_y
          ! zy term
          final(iZeta,ipa,ipb,8) = Two*temp_zy+temp_x
          ! zz term
          final(iZeta,ipa,ipb,9) = Two*temp_zz
        end do

        if (ixb > 0) then
          do iZeta=1,nZeta
            temp_yx = dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3)
            temp_yy = dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),5)
            temp_yz = dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),6)
            temp_zx = -dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2)
            temp_zy = -dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),4)
            temp_zz = -dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),5)

            final(iZeta,ipa,ipb,4) = final(iZeta,ipa,ipb,4)+Two*temp_yx
            final(iZeta,ipa,ipb,5) = final(iZeta,ipa,ipb,5)+Two*temp_yy
            final(iZeta,ipa,ipb,6) = final(iZeta,ipa,ipb,6)+Two*temp_yz
            final(iZeta,ipa,ipb,7) = final(iZeta,ipa,ipb,7)+Two*temp_zx
            final(iZeta,ipa,ipb,8) = final(iZeta,ipa,ipb,8)+Two*temp_zy
            final(iZeta,ipa,ipb,9) = final(iZeta,ipa,ipb,9)+Two*temp_zz
          end do
        end if

        if (iyb > 0) then
          do iZeta=1,nZeta
            temp_xx = -dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),3)
            temp_xy = -dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),5)
            temp_xz = -dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),6)
            temp_zx = dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),1)
            temp_zy = dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),2)
            temp_zz = dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),3)

            final(iZeta,ipa,ipb,1) = final(iZeta,ipa,ipb,1)+Two*temp_xx
            final(iZeta,ipa,ipb,2) = final(iZeta,ipa,ipb,2)+Two*temp_xy
            final(iZeta,ipa,ipb,3) = final(iZeta,ipa,ipb,3)+Two*temp_xz
            final(iZeta,ipa,ipb,7) = final(iZeta,ipa,ipb,7)+Two*temp_zx
            final(iZeta,ipa,ipb,8) = final(iZeta,ipa,ipb,8)+Two*temp_zy
            final(iZeta,ipa,ipb,9) = final(iZeta,ipa,ipb,9)+Two*temp_zz

          end do
        end if

        if (izb > 0) then
          do iZeta=1,nZeta
            temp_xx = dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
            temp_xy = dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),4)
            temp_xz = dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),5)
            temp_yx = -dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1)
            temp_yy = -dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
            temp_yz = -dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3)

            final(iZeta,ipa,ipb,1) = final(iZeta,ipa,ipb,1)+Two*temp_xx
            final(iZeta,ipa,ipb,2) = final(iZeta,ipa,ipb,2)+Two*temp_xy
            final(iZeta,ipa,ipb,3) = final(iZeta,ipa,ipb,3)+Two*temp_xz
            final(iZeta,ipa,ipb,4) = final(iZeta,ipa,ipb,4)+Two*temp_yx
            final(iZeta,ipa,ipb,5) = final(iZeta,ipa,ipb,5)+Two*temp_yy
            final(iZeta,ipa,ipb,6) = final(iZeta,ipa,ipb,6)+Two*temp_yz
          end do
        end if

      end do
    end do

  end do
end do

#ifdef _DEBUGPRINT_
write(6,*) ' In Util3 la,lb=',la,lb
do iElem=1,nElem(la)
  do jElem=1,nElem(lb)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',xx) '
    call RecPrt(Label,' ',final(1,iElem,jElem,1),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',xy) '
    call RecPrt(Label,' ',final(1,iElem,jElem,2),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',xz) '
    call RecPrt(Label,' ',final(1,iElem,jElem,3),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',yx) '
    call RecPrt(Label,' ',final(1,iElem,jElem,4),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',yy) '
    call RecPrt(Label,' ',final(1,iElem,jElem,5),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',yz) '
    call RecPrt(Label,' ',final(1,iElem,jElem,6),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',zx) '
    call RecPrt(Label,' ',final(1,iElem,jElem,7),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',zy) '
    call RecPrt(Label,' ',final(1,iElem,jElem,8),nZeta,1)
    write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,',zz) '
    call RecPrt(Label,' ',final(1,iElem,jElem,9),nZeta,1)
  end do
end do
#endif

return

end subroutine Util3
