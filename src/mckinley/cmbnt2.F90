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

subroutine CmbnT2(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Beta,Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
!***********************************************************************
!                                                                      *
! Object: compute the 2nd derivative of the overlap matrix.            *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl
use Constants, only: One, Two, Four, Six, Half, OneHalf
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: nZeta, la, lb, nHess, IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1), indgrd(0:2,0:1,0:nirrep-1), iu, iv, nOp(2)
real(kind=wp) :: Rnxyz(nZeta,3,0:la+2,0:lb+2), Zeta(nZeta), rKappa(nZeta), rFinal(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), &
                 Alpha(nZeta), Beta(nZeta), Hess(nHess), DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
logical(kind=iwp) :: IfHss(0:1,0:2,0:1,0:2)
integer(kind=iwp) :: ia(3), iax, iay, ib(3), ibx, iby, iCar, ich, iCnt, iCoor, iHess, iIrrep, iMax, ipa, ipb, istab(0:1), iyaMax, &
                     iybMax, iZeta, jCar, jCnt, jCoor, kCoor, nDAO
real(kind=wp) :: Fact, oj, ps
integer(kind=iwp), external :: iPrmt
real(kind=r8), external :: DDot_
! Statement function for Cartesian index
integer(kind=iwp) :: Ind, ixyz, ix, iz, I, itot, i1, i2
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
! Index in the lower triang. local hessian i1 row i2 column
! itot*0  added to avoid compiler warning
I(itot,i1,i2) = itot*0+i1*(i1-1)/2+i2

!iRout = 134
iStab(0) = iu
iStab(1) = iv
!iPrint = nPrint(iRout)
!call GetMem(' Enter CmbnT2','LIST','REAL',iDum,iDum)

do iZeta=1,nZeta
  rKappa(iZeta) = Half*rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)
end do
!if (iPrint >= 99) then
!  call RecPrt(' In CmbnT2: Zeta  ',' ',Zeta,1,nZeta)
!  call RecPrt(' In CmbnT2: rKappa',' ',rKappa,1,nZeta)
!  call RecPrt(' In CmbnT2: Alpha ',' ',Alpha,1,nZeta)
!  call RecPrt(' In CmbnT2: Beta  ',' ',Beta,1,nZeta)
!end if
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
          if (IfHss(0,iCoor-1,0,iCoor-1)) then
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rKappa(iZeta)* &
                                                       ((Two*Alpha(iZeta))**2* &
                                                        ((Two*Beta(iZeta))**2* &
                                                         (Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)+2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                         Six*Beta(iZeta)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))- &
                                                        Two*Alpha(iZeta)* &
                                                        ((Two*Beta(iZeta))**2* &
                                                         (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)+2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                         Six*Beta(iZeta)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))))
              if (lb > 0) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                         ((Two*Alpha(iZeta))**2*Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                          Two*Alpha(iZeta)*(Four*Beta(iZeta))*real(lb,kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(icoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                         ((Two*Alpha(iZeta))**2*real(ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                          Two*Alpha(iZeta)*real(ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(jcoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                         ((Two*Alpha(iZeta))**2*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                          Two*Alpha(iZeta)*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(kcoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                         ((Two*Alpha(iZeta))**2*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)- &
                                                          Two*Alpha(iZeta)*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2))
              end if
              if (ia(iCoor) > 0) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                         Four*Alpha(iZeta)*real(ia(iCoor),kind=wp)* &
                                                         ((Two*Beta(iZeta))**2* &
                                                          (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)+2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                          Six*Beta(iZeta)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                if (lb > 0) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                           Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                           Four*Alpha(iZeta)*real(ia(iCoor),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(iCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                           real(ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                           Four*Alpha(iZeta)*real(ia(iCoor),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                           real(ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                           Four*Alpha(iZeta)*real(ia(iCoor),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                           real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                           Four*Alpha(iZeta)*real(ia(iCoor),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if (ia(iCoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                         real(ia(iCoor)*(ia(iCoor)-1),kind=wp)* &
                                                         ((Two*Beta(iZeta))**2* &
                                                          (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)+2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                          Six*Beta(iZeta)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                if (lb > 0) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                           real(ia(iCoor)*(ia(iCoor)-1),kind=wp)* &
                                                           Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(iCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*(ia(iCoor)-1)*ib(iCoor)*(ib(iCoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)-2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*(ia(iCoor)-1)*ib(jCoor)*(ib(jCoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kCoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*(ia(iCoor)-1)*ib(kCoor)*(ib(kCoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
            end do
          end if
        end do

        ! Integrals like dI/dxdz

        do kCoor=1,3
          iCoor = mod(kCoor,3)+1
          jCoor = mod(iCoor,3)+1
          iMax = max(iCoor,jCoor)
          jCoor = min(iCoor,jCoor)
          iCoor = iMax
          if (IfHss(0,iCoor-1,0,jCoor-1)) then
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rKappa(iZeta)* &
                                                       (Two*Alpha(iZeta))**2* &
                                                       ((Two*Beta(iZeta))**2* &
                                                        (Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)+2)* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)+2)* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                        Six*Beta(iZeta)* &
                                                        Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                        Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              if (lb > 0) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                         (Two*Alpha(iZeta))**2*Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(icoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                         (Two*Alpha(iZeta))**2*real(ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(jcoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                         (Two*Alpha(iZeta))**2*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(kcoor) > 1) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                         (Two*Alpha(iZeta))**2*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
              end if
              if (ia(icoor) > 0) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                         real(ia(icoor),kind=wp)*Two*Alpha(iZeta)* &
                                                         ((Two*Beta(iZeta))**2* &
                                                          (Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)+2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)+2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                          Six*Beta(iZeta)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                if (lb > 0) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(iCoor))*Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(iCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(iCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(iCoor)*ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if (ia(jcoor) > 0) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                         real(ia(jcoor),kind=wp)*Two*Alpha(iZeta)* &
                                                         ((Two*Beta(iZeta))**2* &
                                                          (Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)+2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)+2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                          Six*Beta(iZeta)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                if (lb > 0) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(jCoor))*Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(jCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)*real(ia(jCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           Two*Alpha(iZeta)* &
                                                           real(ia(jCoor),kind=wp)*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if ((ia(iCoor) > 0) .and. (ia(jCoor) > 0)) then
                rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                         real(ia(iCoor)*ia(jCoor),kind=wp)* &
                                                         ((Two*Beta(iZeta))**2* &
                                                          (Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)+2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)+2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))+ &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                          Six*Beta(iZeta)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                if (lb > 0) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                           real(ia(iCoor)*ia(jCoor),kind=wp)*Four*Beta(iZeta)*real(lb,kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*ia(jCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*ia(jCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rFinal(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                           real(ia(iCoor)*ia(jCoor)*ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                           Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                           Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                           Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if

            end do
          end if
        end do
      end do
    end do
  end do
end do

! Trace the Hessian integrals

nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
!if (iPrint >= 99) then
!  call RecPrt(' S(1)',' ',rFinal,nDAO,21)
!  call RecPrt('   D ',' ',DAO,nDAO,1)
!end if
do iIrrep=0,nIrrep-1
  do iCnt=0,1
    do iCar=1,3
      do jCnt=0,1
        do jCar=1,3
          if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) /= 0) then

            ! Accumulate contribution to the Hessian

            ! Get the characteristics of the diff operator

            iCh = ieor(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)

            ! Get the character of the operator in the present irrep

            ps = real(iChTbl(iIrrep,nOp(2))**(iCnt+jCnt),kind=wp)

            ! Get the transf. character of the diff. operator

            ps = ps*real(iPrmt(nOp(2),iCh),kind=wp)

            ! If the over triangular diff. are needed multiply by two instead

            if ((iCnt /= jCnt) .and. (iCar == jCar) .and. (abs(indgrd(iCar-1,iCnt,iIrrep)) == abs(indgrd(jCar-1,jCnt,iIrrep)))) then
              ps = ps*Two
            end if
            iHess = abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
            Fact = real(iStab(iCnt)*iStab(jCnt),kind=wp)/real(nIrrep**2,kind=wp)
            Fact = Fact*ps
            if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) > 0) then
              oj = DDot_(nDAO,DAO,1,rFinal(1,1,1,I(6,iCar,jCar)),1)
              Hess(iHess) = Hess(iHess)+Fact*oj
            else
              Fact = Fact*(-One)**(icnt+jcnt)
              oj = DDot_(nDAO,DAO,1,rFinal(1,1,1,I(6,max(iCar,jCar),min(iCar,jCar))),1)
              Hess(iHess) = Hess(iHess)+Fact*oj
            end if
          end if
        end do
      end do
    end do
  end do
end do

return

end subroutine CmbnT2
