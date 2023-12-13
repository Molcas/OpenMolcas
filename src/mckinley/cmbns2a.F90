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
! Copyright (C) 1994, Anders Bernhardsson                              *
!               1994, Roland Lindh                                     *
!***********************************************************************

subroutine CmbnS2a(Rnxyz,nZeta,la,lb,rKappa,rFinal,Alpha,IfHss,ld)
!***********************************************************************
!                                                                      *
! Object: compute the 2nd derivative of the overlap matrix.            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: C_Ind, iTri, nTri_Elem1
use Constants, only: Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, ld
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+ld,0:lb), rKappa(nZeta), Alpha(nZeta)
real(kind=wp), intent(inout) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6)
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3)
integer(kind=iwp) :: ia(3), iax, iay, ib(3), ibx, iby, iCoor, ii, ij, ipa, ipb, iyaMax, iybMax, jCoor, kCoor

do iax=0,la
  ia(1) = iax
  iyaMax = la-ia(1)
  do ibx=0,lb
    ib(1) = ibx
    iybMax = lb-ib(1)
    do iay=0,iyaMax
      ia(2) = iay
      ia(3) = la-ia(2)-ia(1)
      ipa = C_Ind(la,ia(1),ia(3))
      do iby=0,iybMax
        ib(2) = iby
        ib(3) = lb-ib(2)-ib(1)
        ipb = C_Ind(lb,ib(1),ib(3))

        ! Combine overlap integrals

        ! Integrals like dI/dx1dx1

        do iCoor=1,3
          jCoor = mod(iCoor,3)+1
          kCoor = mod(jCoor,3)+1
          if (IfHss(1,iCoor,1,iCoor)) then
            ii = iTri(iCoor,iCoor)
            rFinal(:,ipa,ipb,ii) = rKappa(:)*((Two*Alpha(:))**2*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))- &
                                              Two*Alpha(:)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (ia(iCoor) > 0) rFinal(:,ipa,ipb,ii) = rFinal(:,ipa,ipb,ii)-rKappa(:)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                                      Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                      Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                      Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ia(iCoor) > 1) rFinal(:,ipa,ipb,ii) = rFinal(:,ipa,ipb,ii)+rKappa(:)*real(ia(iCoor)*(ia(iCoor)-1),kind=wp)* &
                                                      Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                      Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                      Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
          end if
        end do

        ! Integrals like dI/dxdz

        do iCoor=2,3
          do jCoor=1,iCoor-1
            if (IfHss(1,iCoor,1,jCoor)) then
              ij = iTri(iCoor,jCoor)
              rFinal(:,ipa,ipb,ij) = rKappa
              do kCoor=1,3
                if ((kCoor == iCoor) .or. (kCoor == jCoor)) then
                  if (ia(kCoor) > 0) then
                    rFinal(:,ipa,ipb,ij) = rFinal(:,ipa,ipb,ij)*(Two*Alpha(:)*Rnxyz(:,kCoor,ia(kCoor)+1,ib(kCoor))- &
                                                                 real(ia(kCoor),kind=wp)*Rnxyz(:,kCoor,ia(kCoor)-1,ib(kCoor)))
                  else
                    rFinal(:,ipa,ipb,ij) = rFinal(:,ipa,ipb,ij)*Two*Alpha(:)*Rnxyz(:,kCoor,ia(kCoor)+1,ib(kCoor))
                  end if
                else
                  rFinal(:,ipa,ipb,ij) = rFinal(:,ipa,ipb,ij)*Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
                end if
              end do
            end if
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine CmbnS2a
