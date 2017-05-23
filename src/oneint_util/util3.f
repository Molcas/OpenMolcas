************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991,2015, Roland Lindh                                *
*               2015, Lasse Kragh Soerensen                            *
************************************************************************
      SubRoutine Util3(Beta,nZeta,Final,la,lb,Slalbp,Slalb,Slalbm)
************************************************************************
*                                                                      *
* Object: to assemble the orbital angular momentum integrals from the  *
*         derivative integrals dipole integrals.        .              *
*                                                                      *
* Called from: OMQInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Lasse Kragh Soerensen                      *
*             Dept. of Theoretical Chemistry,                          *
*             University of Uppsala, SWEDEN                            *
*             February '15                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
!     Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,9),
!    &        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3,3),
!    &        Slalb (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3),
!    &        Slalbm(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2,3,3),
!    &        Beta(nZeta)
      Real*8  Final (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,9),
     &        Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,6),
     &        Slalb (nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3),
     &        Slalbm(nZeta,(la+1)*(la+2)/2, lb   *(lb+1)/2,6),
     &        Beta(nZeta)
*define _DEBUG_
#ifdef _DEBUG_
      Character*80 Label
#endif
!
! Notice CmbnMP has just 6 components instead of 9!!! (automatically assumes symmetry) Well fuck you CmbnMP
! This means Slalbp and Slalbm in reality only has 6 components....
! XX = 1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5 and ZZ=6
! There are only six components since zy d/di = yz d/di
! We still keeps the 9 components in final
!
*
*     Statement function for cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
#ifdef _DEBUG_
      nElem(ix) = (ix+1)*(ix+2)/2
#endif
*
#ifdef _DEBUG_
       Write (6,*) ' In Util3 la,lb=',la,lb
       Call RecPrt('Beta',' ',Beta,nZeta,1)
       Do ia = 1, nElem(la)
          Do ib = 1, nElem(lb+1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'xx)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,1),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'xy)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'xz)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'yx)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,2),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'yy)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,4),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'yz)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,5),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'zx)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,3),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'zy)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,5),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalbp(',ia,',',ib,'zz)'
             Call RecPrt(Label,' ',Slalbp(1,ia,ib,6),nZeta,1)
          End Do
       End Do
       Do ia = 1, nElem(la)
          Do ib = 1, nElem(lb)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalb (',ia,',',ib,'x)'
             Call RecPrt(Label,' ',Slalb (1,ia,ib,1),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalb (',ia,',',ib,'y)'
             Call RecPrt(Label,' ',Slalb (1,ia,ib,2),nZeta,1)
             Write (Label,'(A,I2,A,I2,A)')
     &              ' Slalb (',ia,',',ib,'z)'
             Call RecPrt(Label,' ',Slalb (1,ia,ib,3),nZeta,1)
          End Do
       End Do
       If (lb.gt.0) Then
          Do ia = 1, nElem(la)
             Do ib = 1, nElem(lb-1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'xx)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'xy)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'xz)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'yx)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'yy)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,4),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'yz)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'zx)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'zy)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Slalbm(',ia,',',ib,'zz)'
                Call RecPrt(Label,' ',Slalbm(1,ia,ib,6),nZeta,1)
             End Do
          End Do
       End If
#endif
*
      Do 10 ixa = la, 0, -1
         Do 11 iya = la-ixa, 0, -1
            iza = la-ixa-iya
            ipa = Ind(la,ixa,iza)
*
      Do 20 ixb = lb, 0, -1
         Do 21 iyb = lb-ixb, 0, -1
            izb = lb-ixb-iyb
            ipb = Ind(lb,ixb,izb)
*
            Do 30 iZeta = 1, nZeta
!              temp_xx = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1,2)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),1,3) )
!              temp_yx = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2,2)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),2,3) )
!              temp_zx = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3,2)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),3,3) )
!              temp_xy = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1,3)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1,1) )
!              temp_yy = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2,3)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2,1) )
!              temp_zy = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3,3)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3,1) )
!              temp_xz = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),1,1)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1,2) )
!              temp_yz = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),2,1)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2,2) )
!              temp_zz = Two*Beta(iZeta) * (
!    &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),3,1)
!    &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3,2) )
               temp_xx = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),3) )
               temp_xy = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),4)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),5) )
               temp_xz = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),5)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),6) )
               temp_yx = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1) )
               temp_yy = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),5)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2) )
               temp_yz = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),6)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3) )
               temp_zx = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),1)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2) )
               temp_zy = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),2)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),4) )
               temp_zz = Two*Beta(iZeta) * (
     &                  Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb  ),3)
     &                 -Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),5) )
*
               temp_x = Slalb (iZeta,ipa,Ind(lb,ixb,izb),1)
               temp_y = Slalb (iZeta,ipa,Ind(lb,ixb,izb),2)
               temp_z = Slalb (iZeta,ipa,Ind(lb,ixb,izb),3)
*
*              xx term
               Final(iZeta,ipa,ipb,1) = Two * temp_xx
*              xy term
               Final(iZeta,ipa,ipb,2) = Two * temp_xy - temp_z
*              xz term
               Final(iZeta,ipa,ipb,3) = Two * temp_xz + temp_y
*              yx term
               Final(iZeta,ipa,ipb,4) = Two * temp_yx + temp_z
*              yy term
               Final(iZeta,ipa,ipb,5) = Two * temp_yy
*              yz term
               Final(iZeta,ipa,ipb,6) = Two * temp_yz - temp_x
*              zx term
               Final(iZeta,ipa,ipb,7) = Two * temp_zx - temp_y
*              zy term
               Final(iZeta,ipa,ipb,8) = Two * temp_zy + temp_x
*              zz term
               Final(iZeta,ipa,ipb,9) = Two * temp_zz
 30         Continue
*
            If (ixb.gt.0) Then
               Do 31 iZeta = 1, nZeta
!              temp_xy = Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1,3)
!              temp_yy = Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2,3)
!              temp_zy = Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3,3)
!              temp_xz =-Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1,2)
!              temp_yz =-Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2,2)
!              temp_zz =-Dble(ixb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3,2)
               temp_yx = Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3)
               temp_yy = Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),5)
               temp_yz = Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),6)
               temp_zx =-Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2)
               temp_zy =-Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),4)
               temp_zz =-Dble(ixb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),5)
*
               Final(iZeta,ipa,ipb,4) = Final(iZeta,ipa,ipb,4)
     &                                                   + Two * temp_yx
               Final(iZeta,ipa,ipb,5) = Final(iZeta,ipa,ipb,5)
     &                                                   + Two * temp_yy
               Final(iZeta,ipa,ipb,6) = Final(iZeta,ipa,ipb,6)
     &                                                   + Two * temp_yz
               Final(iZeta,ipa,ipb,7) = Final(iZeta,ipa,ipb,7)
     &                                                   + Two * temp_zx
               Final(iZeta,ipa,ipb,8) = Final(iZeta,ipa,ipb,8)
     &                                                   + Two * temp_zy
               Final(iZeta,ipa,ipb,9) = Final(iZeta,ipa,ipb,9)
     &                                                   + Two * temp_zz
 31            Continue
            End If
*
            If (iyb.gt.0) Then
               Do 32 iZeta = 1, nZeta
!              temp_xx =-Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),1,3)
!              temp_yx =-Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),2,3)
!              temp_zx =-Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),3,3)
!              temp_xz = Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),1,1)
!              temp_yz = Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),2,1)
!              temp_zz = Dble(iyb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),3,1)
               temp_xx =-Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),3)
               temp_xy =-Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),5)
               temp_xz =-Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),6)
               temp_zx = Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),1)
               temp_zy = Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),2)
               temp_zz = Dble(iyb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb  ),3)
*
               Final(iZeta,ipa,ipb,1) = Final(iZeta,ipa,ipb,1)
     &                                                   + Two * temp_xx
               Final(iZeta,ipa,ipb,2) = Final(iZeta,ipa,ipb,2)
     &                                                   + Two * temp_xy
               Final(iZeta,ipa,ipb,3) = Final(iZeta,ipa,ipb,3)
     &                                                   + Two * temp_xz
               Final(iZeta,ipa,ipb,7) = Final(iZeta,ipa,ipb,7)
     &                                                   + Two * temp_zx
               Final(iZeta,ipa,ipb,8) = Final(iZeta,ipa,ipb,8)
     &                                                   + Two * temp_zy
               Final(iZeta,ipa,ipb,9) = Final(iZeta,ipa,ipb,9)
     &                                                   + Two * temp_zz
*
 32            Continue
            End If
*
            If (izb.gt.0) Then
               Do 33 iZeta = 1, nZeta
!              temp_xx = Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1,2)
!              temp_yx = Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2,2)
!              temp_zx = Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3,2)
!              temp_xy =-Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1,1)
!              temp_yy =-Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2,1)
!              temp_zy =-Dble(izb) *
!    &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3,1)
               temp_xx = Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
               temp_xy = Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),4)
               temp_xz = Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),5)
               temp_yx =-Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1)
               temp_yy =-Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
               temp_yz =-Dble(izb) *
     &                  Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3)
*
               Final(iZeta,ipa,ipb,1) = Final(iZeta,ipa,ipb,1)
     &                                                   + Two * temp_xx
               Final(iZeta,ipa,ipb,2) = Final(iZeta,ipa,ipb,2)
     &                                                   + Two * temp_xy
               Final(iZeta,ipa,ipb,3) = Final(iZeta,ipa,ipb,3)
     &                                                   + Two * temp_xz
               Final(iZeta,ipa,ipb,4) = Final(iZeta,ipa,ipb,4)
     &                                                   + Two * temp_yx
               Final(iZeta,ipa,ipb,5) = Final(iZeta,ipa,ipb,5)
     &                                                   + Two * temp_yy
               Final(iZeta,ipa,ipb,6) = Final(iZeta,ipa,ipb,6)
     &                                                   + Two * temp_yz
 33            Continue
            End If
*
 21      Continue
 20   Continue
*
 11      Continue
 10   Continue
*
#ifdef _DEBUG_
          Write (6,*) ' In Util3 la,lb=',la,lb
          Do 300 iElem = 1, nElem(la)
             Do 310 jElem = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xx) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xy) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',xz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',yx) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,4),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',yy) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',yz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,6),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',zx) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,7),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',zy) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,8),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Final (',iElem,',',jElem,',zz) '
                Call RecPrt(Label,' ',Final(1,iElem,jElem,9),nZeta,1)
 310         Continue
 300      Continue
#endif
*
      Return
      End
