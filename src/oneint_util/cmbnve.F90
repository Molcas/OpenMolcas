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

subroutine CmbnVe(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,final,nComp,Vxyz)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), rKappa(nZeta), Rnxyz(nZeta,3,0:la,0:lb+1,0:lr), &
       Vxyz(nZeta,3,0:la,0:lb)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

!iRout = 161
!iPrint = nPrint(iRout)

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
        !if (iPrint >= 99) then
        !  write(6,*) ixa,iya,iza,ixb,iyb,izb
        !  write(6,*) ipa,ipb
        !end if

        ! Combine integrals

        do iZeta=1,nZeta
          Fact = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)
          final(iZeta,ipa,ipb,1) = Fact*Vxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb,0)*Rnxyz(iZeta,3,iza,izb,0)
          final(iZeta,ipa,ipb,2) = Fact*Rnxyz(iZeta,1,ixa,ixb,0)*Vxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb,0)
          final(iZeta,ipa,ipb,3) = Fact*Rnxyz(iZeta,1,ixa,ixb,0)*Rnxyz(iZeta,2,iya,iyb,0)*Vxyz(iZeta,3,iza,izb)
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnVe
