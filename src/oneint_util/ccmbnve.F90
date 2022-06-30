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

subroutine CCmbnVe(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,nComp,Vxyz,KVector,P)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January  91                                              *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Half, Quart, OneHalf, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp
complex(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1), Vxyz(nZeta,3,0:la,0:lb,2)
real(kind=wp), intent(in) :: Zeta(nZeta), rKappa(nZeta), kVector(3), P(nZeta,3)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp)
#include "print.fh"
integer(kind=iwp) :: ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta
real(kind=wp) :: Fact, k_dot_P, kModQ, rTemp
complex(kind=wp) :: Temp1, Temp2, Tempm, Tempp

iRout = 161
iPrint = nPrint(iRout)

kModQ = Quart*(kVector(1)**2+kVector(2)**2+kVector(3)**2)

do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)

        ! Combine integrals

        do iZeta=1,nZeta

          ! Put in the correct prefactors

          rTemp = kModQ/Zeta(iZeta)
          Fact = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)*exp(-rTemp)
          k_dot_P = kVector(1)*P(iZeta,1)+kVector(2)*P(iZeta,2)+kVector(3)*P(iZeta,3)
          Temp1 = Vxyz(iZeta,1,ixa,ixb,1)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
          Temp2 = Vxyz(iZeta,1,ixa,ixb,2)*Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
          Tempp = exp(Onei*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(Onei*k_dot_P)*Fact*(Temp1-Temp2)*Half
          rFinal(iZeta,ipa,ipb,1) = real(Tempp)
          rFinal(iZeta,ipa,ipb,4) = real(Tempm)
          rFinal(iZeta,ipa,ipb,7) = aimag(Tempp)
          rFinal(iZeta,ipa,ipb,10) = aimag(Tempm)
          Temp1 = Rnxyz(iZeta,1,ixa,ixb)*Vxyz(iZeta,2,iya,iyb,1)*Rnxyz(iZeta,3,iza,izb)
          Temp2 = Rnxyz(iZeta,1,ixa,ixb)*Vxyz(iZeta,2,iya,iyb,2)*Rnxyz(iZeta,3,iza,izb)
          Tempp = exp(Onei*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(Onei*k_dot_P)*Fact*(Temp1-Temp2)*Half
          rFinal(iZeta,ipa,ipb,2) = real(Tempp)
          rFinal(iZeta,ipa,ipb,5) = real(Tempm)
          rFinal(iZeta,ipa,ipb,8) = aimag(Tempp)
          rFinal(iZeta,ipa,ipb,11) = aimag(Tempm)
          Temp1 = Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Vxyz(iZeta,3,iza,izb,1)
          Temp2 = Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*Vxyz(iZeta,3,iza,izb,2)
          Tempp = exp(Onei*k_dot_P)*Fact*(Temp1+Temp2)*Half
          Tempm = exp(Onei*k_dot_P)*Fact*(Temp1-Temp2)*Half
          rFinal(iZeta,ipa,ipb,3) = real(Tempp)
          rFinal(iZeta,ipa,ipb,6) = real(Tempm)
          rFinal(iZeta,ipa,ipb,9) = aimag(Tempp)
          rFinal(iZeta,ipa,ipb,12) = aimag(Tempm)
        end do
        if (iPrint >= 99) then
          write(u6,*) '(',ixa,iya,iza,ixb,iyb,izb,')'
          write(u6,*) 'x-component'
          write(u6,*) rFinal(1,ipa,ipb,1)
          write(u6,*) rFinal(1,ipa,ipb,4)
          write(u6,*) rFinal(1,ipa,ipb,7)
          write(u6,*) rFinal(1,ipa,ipb,10)
          write(u6,*) 'y-component'
          write(u6,*) rFinal(1,ipa,ipb,2)
          write(u6,*) rFinal(1,ipa,ipb,5)
          write(u6,*) rFinal(1,ipa,ipb,8)
          write(u6,*) rFinal(1,ipa,ipb,11)
          write(u6,*) 'z-component'
          write(u6,*) rFinal(1,ipa,ipb,3)
          write(u6,*) rFinal(1,ipa,ipb,6)
          write(u6,*) rFinal(1,ipa,ipb,9)
          write(u6,*) rFinal(1,ipa,ipb,12)
        end if

      end do
    end do
  end do
end do

return

end subroutine CCmbnVe
