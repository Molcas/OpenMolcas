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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnS1_mck(Rnxyz,nZeta,la,lb,Zeta,rKappa,final,Alpha,Beta,IfGrad,nOp)
!***********************************************************************
!                                                                      *
! Object: compute the gradient of the overlap matrix.                  *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!             October '91.                                             *
!             Anders Bernhardsson                                      *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!              95.                                                     *
!***********************************************************************

implicit real*8(A-H,O-Z)
!#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,1), Zeta(nZeta), rKappa(nZeta), Beta(nZeta), Rnxyz(nZeta,3,0:la+1,0:lb+1), &
       Alpha(nZeta)
logical IfGrad(3,2)
integer nOp(2)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

!iRout = 134
!iPrint = nPrint(iRout)
!call GetMem(' Enter CmbnS1_mck','LIST','REAL',iDum,iDum)

!ii = la*(la+1)*(la+2)/6
!jj = lb*(lb+1)*(lb+2)/6
exp32 = -Three/Two
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**exp32
end do
!if (iPrint >= 99) then
!  call RecPrt(' In CmbnS1_mck: Zeta  ',' ',Zeta,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: rKappa',' ',rKappa,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: Alpha ',' ',Alpha,1,nZeta)
!  call RecPrt(' In CmbnS1_mck: Beta  ',' ',Beta,1,nZeta)
!end if
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

        ! Combine overlap integrals

        tTwo = Two
        !write(6,*) ' papb=', papb
        if (IfGrad(1,1)) then
          if (ixa > 0) then
            xa = dble(-ixa)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*(tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)+xa*Rnxyz(iZeta,1,ixa-1,ixb))* &
                                       Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb)*Rnxyz(iZeta,2,iya,iyb)* &
                                       Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(1,2)) then
          if (ixb > 0) then
            xb = dble(-ixb)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)+xb*Rnxyz(iZeta,1,ixa,ixb-1))* &
                                       Rnxyz(iZeta,2,iya,iyb)*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1)*Rnxyz(iZeta,2,iya,iyb)* &
                                       Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,1)) then
          if (iya > 0) then
            ya = dble(-iya)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*(tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)+ &
                                       ya*Rnxyz(iZeta,2,iya-1,iyb))*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb)* &
                                       Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(2,2)) then
          if (iyb > 0) then
            yb = dble(-iyb)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)+ &
                                       yb*Rnxyz(iZeta,2,iya,iyb-1))*Rnxyz(iZeta,3,iza,izb)
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1)* &
                                       Rnxyz(iZeta,3,iza,izb)
            end do
          end if
        end if
        if (IfGrad(3,1)) then
          if (iza > 0) then
            za = dble(-iza)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(tTwo*Alpha(iZeta)* &
                                       Rnxyz(iZeta,3,iza+1,izb)+za*Rnxyz(iZeta,3,iza-1,izb))
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*tTwo*Alpha(iZeta)* &
                                       Rnxyz(iZeta,3,iza+1,izb)
            end do
          end if
        end if
        if (IfGrad(3,2)) then
          if (izb > 0) then
            zb = dble(-izb)
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*(tTwo*Beta(iZeta)* &
                                       Rnxyz(iZeta,3,iza,izb+1)+zb*Rnxyz(iZeta,3,iza,izb-1))
            end do
          else
            do iZeta=1,nZeta
              !final(iZeta,ipa,ipb,1) = papb*rKappa(iZeta)*
              final(iZeta,ipa,ipb,1) = rKappa(iZeta)*Rnxyz(iZeta,1,ixa,ixb)*Rnxyz(iZeta,2,iya,iyb)*tTwo*Beta(iZeta)* &
                                       Rnxyz(iZeta,3,iza,izb+1)
            end do
          end if
        end if

      end do
    end do
  end do
end do

!call GetMem(' Exit CmbnS1_mck','LIST','REAL',iDum,iDum)
return

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(nOp)
end if

end subroutine CmbnS1_mck
