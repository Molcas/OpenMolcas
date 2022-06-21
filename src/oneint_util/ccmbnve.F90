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

subroutine CCmbnVe(Rnxyz,nZeta,la,lb,Zeta,rKappa,final,nComp,Vxyz,KVector,P)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January  91                                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), rKappa(nZeta), rTemp, Fact, KVector(3), P(nZeta,3), k_dot_P
complex*16 Rnxyz(nZeta,3,0:la+1,0:lb+1), Vxyz(nZeta,3,0:la,0:lb,2), Temp1, Temp2, Tempp, Tempm, i
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

iRout = 161
iPrint = nPrint(iRout)

i = (0.0d0,1.0d0)
do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        ! Combine integrals

        do iZeta=1,nZeta

          ! Put in the correct prefactors

          rTemp = KVector(1)**2+kVector(2)**2+kVector(3)**2
          rTemp = rTemp/(Four*Zeta(iZeta))
          Fact = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)*exp(-rTemp)
          k_dot_P = kVector(1)*P(iZeta,1)+kVector(2)*P(iZeta,2)+kVector(3)*P(iZeta,3)
          Temp1 = Vxyz(iZeta,1,ixa,ixb,1)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
          Temp2 = Vxyz(iZeta,1,ixa,ixb,2)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
          Tempp = exp(i*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(i*k_dot_P)*Fact*(Temp1-Temp2)*Half
          final(iZeta,ipa,ipb,1) = dble(Tempp)
          final(iZeta,ipa,ipb,4) = dble(Tempm)
          final(iZeta,ipa,ipb,7) = DIMAG(Tempp)
          final(iZeta,ipa,ipb,10) = DIMAG(Tempm)
          Temp1 = Rnxyz(iZeta,1,ixa,ixb)*Vxyz(iZeta,2,iya,iyb,1)*Rnxyz(iZeta,3,iza,izb)
          Temp2 = Rnxyz(iZeta,1,ixa,ixb)*Vxyz(iZeta,2,iya,iyb,2)*Rnxyz(iZeta,3,iza,izb)
          Tempp = exp(i*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(i*k_dot_P)*Fact*(Temp1-Temp2)*Half
          final(iZeta,ipa,ipb,2) = dble(Tempp)
          final(iZeta,ipa,ipb,5) = dble(Tempm)
          final(iZeta,ipa,ipb,8) = DIMAG(Tempp)
          final(iZeta,ipa,ipb,11) = DIMAG(Tempm)
          Temp1 = Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Vxyz(iZeta,3,iza,izb,1)
          Temp2 = Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Vxyz(iZeta,3,iza,izb,2)
          Tempp = exp(i*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(i*k_dot_P)*Fact*(Temp1-Temp2)*Half
          final(iZeta,ipa,ipb,3) = dble(Tempp)
          final(iZeta,ipa,ipb,6) = dble(Tempm)
          final(iZeta,ipa,ipb,9) = DIMAG(Tempp)
          final(iZeta,ipa,ipb,12) = DIMAG(Tempm)
        end do
        if (iPrint >= 99) then
          write(6,*) '(',ixa,iya,iza,ixb,iyb,izb,')'
          write(6,*) 'x-component'
          write(6,*) final(1,ipa,ipb,1)
          write(6,*) final(1,ipa,ipb,4)
          write(6,*) final(1,ipa,ipb,7)
          write(6,*) final(1,ipa,ipb,10)
          write(6,*) 'y-component'
          write(6,*) final(1,ipa,ipb,2)
          write(6,*) final(1,ipa,ipb,5)
          write(6,*) final(1,ipa,ipb,8)
          write(6,*) final(1,ipa,ipb,11)
          write(6,*) 'z-component'
          write(6,*) final(1,ipa,ipb,3)
          write(6,*) final(1,ipa,ipb,6)
          write(6,*) final(1,ipa,ipb,9)
          write(6,*) final(1,ipa,ipb,12)
        end if

      end do
    end do
  end do
end do

return

end subroutine CCmbnVe
