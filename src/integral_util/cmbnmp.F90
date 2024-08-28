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

subroutine CmbnMP(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,final,nComp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

implicit none
integer nZeta, la, lb, nComp, lr
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Zeta(nZeta), rKappa(nZeta), Rnxyz(nZeta,3,0:la,0:lb,0:lr)
integer ixa, ixb, iya, iyb, iza, izb, iyaMax, iybMax, ipa, ipb, ix, iy, iz
integer ixyz, iComp, iZeta, Ind
real*8 Fact
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

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
        !  write(u6,*) ixa,iya,iza,ixb,iyb,izb
        !  write(u6,*) ipa,ipb
        !end if

        ! Combine multipole moment integrals

        iComp = 0
        do ix=lr,0,-1
          do iy=lr-ix,0,-1
            iz = lr-ix-iy
            iComp = iComp+1
            !write(u6,*) ix,iy,iz,iComp
            do iZeta=1,nZeta
              Fact = rKappa(iZeta)*1/sqrt(Zeta(iZeta)**3)
              !Fact = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)
              final(iZeta,ipa,ipb,iComp) = Fact*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)*Rnxyz(iZeta,3,iza,izb,iz)
            end do
          end do
        end do

      end do
    end do
  end do
end do

return

end subroutine CmbnMP
