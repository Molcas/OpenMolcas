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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine TraPAB(nZeta,la,lb,AB,GInt,jSum,rKappa,Fac1,Fac2,Fac3,Fac4,Fac5,A,B,P)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Index_Functions, only: C_Ind3, C3_Ind3, nTri_Elem1
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri3_Elem1
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, jSum
real(kind=wp), intent(out) :: AB(nZeta,nTri_Elem1(la),nTri_Elem1(lb)), Fac1(nZeta), Fac2(nZeta), Fac3(nZeta), Fac4(nZeta), &
                              Fac5(nZeta)
real(kind=wp), intent(inout) :: GInt(nZeta,jSum)
real(kind=wp), intent(in) :: rKappa(nZeta), A(3), B(3), P(nZeta,3)
integer(kind=iwp) :: i, ia, iax, iay, iaz, ibx, iby, ibz, igx, igy, igz, ipa, ipb, ipg, iTrgt, ix, ixa, ixas, ixb, ixbs, iy, iya, &
                     iyaMax, iyas, iyb, iybMax, iybs, iz, iza, izas, izb, izbs, iZeta, jx, jy, jz, kOff, lOff
real(kind=wp) :: Ax, Ay, Az, Bx, By, Bz

#ifdef _DEBUGPRINT_
call RecPrt(' In TraPAB: GInt',' ',GInt,nZeta,jSum)
call RecPrt(' In TraPAB: P   ',' ',P,nZeta,3)
#endif

! Initialize

AB(:,:,:) = Zero

! Remove redundant elements in GInt. This is done in place.

kOff = 4
do i=2,la+lb

  do ix=i,0,-1
    do iy=i-ix,0,-1
      iz = i-ix-iy
      jx = ix
      jy = iy
      jz = iz

      lOff = 0
      do ia=1,i-1
        if (jz /= 0) then
          lOff = 3*(2+lOff)
          jz = jz-1
        else if (jy /= 0) then
          lOff = 3*(1+lOff)
          jy = jy-1
        else
          lOff = 3*lOff
          jx = jx-1
        end if
      end do
      if (jz == 1) lOff = lOff+3
      if (jy == 1) lOff = lOff+2
      if (jx == 1) lOff = lOff+1

      iTrgt = C3_Ind3(ix,iy,iz)
      !write(u6,*) ' ix,iy,iz,kOff,lOff,iTrgt=',ix,iy,iz,kOff,lOff,iTrgt
      GInt(:,iTrgt) = GInt(:,kOff+lOff)

    end do
  end do

  kOff = kOff+3**i
end do
#ifdef _DEBUGPRINT_
call RecPrt(' In TraPAB: GInt(unique)',' ',GInt,nZeta,nTri3_Elem1(la+lb))
#endif

! Loop over the elements of the basis functions on A and B

do ixa=la,0,-1
  iyaMax = la-ixa
  do iya=iyaMax,0,-1
    iza = la-ixa-iya
    ipa = C_Ind3(ixa,iya,iza)
    !write(u6,*) ' ipa,ixa,iya,iza=',ipa,ixa,iya,iza

    do ixb=lb,0,-1
      iybMax = lb-ixb
      do iyb=iybMax,0,-1
        izb = lb-ixb-iyb
        ipb = C_Ind3(ixb,iyb,izb)
        !write(u6,*) ' ipb,ixb,iyb,izb=',ipb,ixb,iyb,izb

        ! Loop over the elements of functions at P

        do ixas=0,ixa
          call Binom(ixa,ixas,iAx)
          Ax = real(iAx,kind=wp)
          do iZeta=1,nZeta
            if ((ixa-ixas) == 0) then
              Fac1(iZeta) = rKappa(iZeta)*Ax
            else
              Fac1(iZeta) = rKappa(iZeta)*Ax*(P(iZeta,1)-A(1))**(ixa-ixas)
            end if
          end do
          do iyas=0,iya
            call Binom(iya,iyas,iAy)
            Ay = real(iAy,kind=wp)
            do iZeta=1,nZeta
              if ((iya-iyas) == 0) then
                Fac2(iZeta) = Fac1(iZeta)*Ay
              else
                Fac2(iZeta) = Fac1(iZeta)*Ay*(P(iZeta,2)-A(2))**(iya-iyas)
              end if
            end do
            do izas=0,iza
              call Binom(iza,izas,iAz)
              Az = real(iAz,kind=wp)
              do iZeta=1,nZeta
                if ((iza-izas) == 0) then
                  Fac3(iZeta) = Fac2(iZeta)*Az
                else
                  Fac3(iZeta) = Fac2(iZeta)*Az*(P(iZeta,3)-A(3))**(iza-izas)
                end if
              end do

              do ixbs=0,ixb
                call Binom(ixb,ixbs,iBx)
                Bx = real(iBx,kind=wp)
                do iZeta=1,nZeta
                  if ((ixb-ixbs) == 0) then
                    Fac4(iZeta) = Fac3(iZeta)*Bx
                  else
                    Fac4(iZeta) = Fac3(iZeta)*Bx*(P(iZeta,1)-B(1))**(ixb-ixbs)
                  end if
                end do
                igx = ixas+ixbs
                do iybs=0,iyb
                  call Binom(iyb,iybs,iBy)
                  By = real(iBy,kind=wp)
                  do iZeta=1,nZeta
                    if ((iyb-iybs) == 0) then
                      Fac5(iZeta) = Fac4(iZeta)*By
                    else
                      Fac5(iZeta) = Fac4(iZeta)*By*(P(iZeta,2)-B(2))**(iyb-iybs)
                    end if
                  end do
                  igy = iyas+iybs
                  do izbs=0,izb
                    call Binom(izb,izbs,iBz)
                    Bz = real(iBz,kind=wp)
                    igz = izas+izbs
                    ipg = C3_Ind3(igx,igy,igz)
                    !write(u6,*) ' ipg,igx,igy,igz=', ipg,igx,igy,igz

                    do iZeta=1,nZeta
                      if ((izb-izbs) == 0) then
                        AB(iZeta,ipa,ipb) = AB(iZeta,ipa,ipb)+Fac5(iZeta)*GInt(iZeta,ipg)*Bz
                      else
                        AB(iZeta,ipa,ipb) = AB(iZeta,ipa,ipb)+Fac5(iZeta)*(P(iZeta,3)-B(3))**(izb-izbs)*GInt(iZeta,ipg)*Bz
                      end if
                    end do

                  end do
                end do
              end do

            end do
          end do
        end do

      end do
    end do

  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In TraPAB: AB',' ',AB,nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
#endif

end subroutine TraPAB
