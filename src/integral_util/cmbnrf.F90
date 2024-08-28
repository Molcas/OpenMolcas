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
! Copyright (C) 1991,1992, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CmbnRF(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,final,nComp,Fact,Temp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!***********************************************************************

use Constants, only: Two, Three
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer nZeta, la, lb, nComp, lr
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), rKappa(nZeta), Fact(nZeta), Temp(nZeta), &
       Rnxyz(nZeta,3,0:la,0:lb,0:lr)
#ifdef _DEBUGPRINT_
integer :: nFinal
#endif
integer ixa, ixb, iya, iyb, iza, izb, ipa, ipb, iyaMax, iybMax, iZeta, iy, ir, iComp
integer ixyz, ix, iz, Ind, iOff
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

do iZeta=1,nZeta
  Fact(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)
end do
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
#       ifdef _DEBUGPRINT_
        write(u6,*) ixa,iya,iza,ixb,iyb,izb
        write(u6,*) ipa,ipb
#       endif

        ! Combine multipole moment integrals

        do ix=0,lr
          do iy=0,lr-ix
            do iZeta=1,nZeta
              Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
            end do
            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              do iZeta=1,nZeta
                final(iZeta,ipa,ipb,iComp) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
              end do
            end do
          end do
        end do

      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
nFinal = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
call RecPrt('Final',' ',final,nFinal,nComp)
#endif

return

end subroutine CmbnRF
