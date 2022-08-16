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

subroutine CmbnS2a(Rnxyz,nZeta,la,lb,rKappa,final,Alpha,IfHss,ld)
!***********************************************************************
!                                                                      *
! Object: compute the 2nd derivative  of the overlap matrix.           *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
!#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), rKappa(nZeta), Rnxyz(nZeta,3,0:la+ld,0:lb), Alpha(nZeta)
logical IfHss(4,3,4,3)
integer ia(3), ib(3)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
! Index in the triang. local hessian
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

do iax=0,la
  ia(1) = iax
  iyaMax = la-ia(1)
  do ibx=0,lb
    ib(1) = ibx
    iybMax = lb-ib(1)
    do iay=0,iyaMax
      ia(2) = iay
      ia(3) = la-ia(2)-ia(1)
      ipa = Ind(la,ia(1),ia(3))
      do iby=0,iybMax
        ib(2) = iby
        ib(3) = lb-ib(2)-ib(1)
        ipb = Ind(lb,ib(1),ib(3))

        ! Combine overlap integrals

        ! Integrals like dI/dx1dx1

        do iCoor=1,3
          jCoor = mod(iCoor,3)+1
          kCoor = mod(jCoor,3)+1
          if (IfHss(1,iCoor,1,iCoor)) then
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,itri(iCoor,iCoor)) = rKappa(iZeta)*((Two*Alpha(iZeta))**2* &
                                                                      Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                      Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                      Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                                      Two*Alpha(iZeta)*Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                                      Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                      Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              if (ia(iCoor) > 0) then
                final(iZeta,ipa,ipb,itri(iCoor,iCoor)) = final(iZeta,ipa,ipb,itri(iCoor,iCoor))- &
                                                         rKappa(iZeta)*(Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                                        Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                                        Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ia(iCoor) > 1) then
                final(iZeta,ipa,ipb,itri(iCoor,iCoor)) = final(iZeta,ipa,ipb,itri(iCoor,iCoor))+ &
                                                         rKappa(iZeta)*(dble(ia(iCoor)*(ia(iCoor)-1))* &
                                                                        Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                                        Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
            end do
          end if
        end do

        ! Integrals like dI/dxdz

        do iCoor=2,3
          do jCoor=1,iCoor-1
            if (IfHss(1,iCoor,1,jCoor)) then
              do kCoor=1,3
                do iZeta=1,nZeta
                  if (kCoor == 1) then
                    final(iZeta,ipa,ipb,itri(iCoor,jCoor)) = rKappa(iZeta)
                  end if
                  if ((kCoor == iCoor) .or. (kCoor == jCoor)) then
                    rIc = Two*Alpha(iZeta)*Rnxyz(iZeta,kCoor,ia(kCoor)+1,ib(kCoor))

                    if (ia(kCoor) > 0) rIc = rIc-dble(ia(kCoor))*Rnxyz(iZeta,kCoor,ia(kCoor)-1,ib(kCoor))

                    final(iZeta,ipa,ipb,itri(iCoor,jCoor)) = final(iZeta,ipa,ipb,itri(iCoor,jCoor))*rIc
                  else
                    final(iZeta,ipa,ipb,itri(iCoor,jCoor)) = final(iZeta,ipa,ipb,itri(iCoor,jCoor))* &
                                                             Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                  end if
                end do
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
