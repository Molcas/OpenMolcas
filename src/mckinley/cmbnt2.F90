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

subroutine CmbnT2(Rnxyz,nZeta,la,lb,Zeta,rKappa,final,Alpha,Beta,Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
!***********************************************************************
!                                                                      *
! Object: compute the 2nd derivative  of the overlap matrix.           *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6), Zeta(nZeta), rKappa(nZeta), Beta(nZeta), Rnxyz(nZeta,3,0:la+2,0:lb+2), &
       Alpha(nZeta), Hess(nHess), DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
logical IfHss(0:1,0:2,0:1,0:2)
integer IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1), istab(0:1), nOp(2), ia(3), ib(3), indgrd(0:2,0:1,0:nirrep-1)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
! Index in the lower triang. local hessian i1 row i2 column
! itot*0  added to avoid compiler warning
I(itot,i1,i2) = itot*0+i1*(i1-1)/2+i2

!iRout = 134
iStab(0) = iu
iStab(1) = iv
!iPrint = nPrint(iRout)
!call GetMem(' Enter CmbnT2','LIST','REAL',iDum,iDum)

exp32 = -Three/Two
do iZeta=1,nZeta
  rKappa(iZeta) = half*rKappa(iZeta)*Zeta(iZeta)**exp32
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
              final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = rKappa(iZeta)* &
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
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                        ((Two*Alpha(iZeta))**2*Four*Beta(iZeta)*dble(lb)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                         Two*Alpha(iZeta)*(Four*Beta(iZeta))*dble(lb)* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(icoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                        ((Two*Alpha(iZeta))**2*dble(ib(icoor)*(ib(icoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)-2)* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                         Two*Alpha(iZeta)*dble(ib(icoor)*(ib(icoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(jcoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                        ((Two*Alpha(iZeta))**2*dble(ib(jcoor)*(ib(jcoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                         Two*Alpha(iZeta)*dble(ib(jcoor)*(ib(jcoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ib(kcoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                        ((Two*Alpha(iZeta))**2*dble(ib(kcoor)*(ib(kcoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)- &
                                                         Two*Alpha(iZeta)*dble(ib(kcoor)*(ib(kcoor)-1))* &
                                                         Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                         Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                         Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2))
              end if
              if (ia(iCoor) > 0) then
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                        Four*Alpha(iZeta)*dble(ia(iCoor))* &
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
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                          Four*Beta(iZeta)*dble(lb)*Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(iCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                          dble(ib(icoor)*(ib(icoor)-1))*Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                          dble(ib(jcoor)*(ib(jcoor)-1))*Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                          dble(ib(kcoor)*(ib(kcoor)-1))*Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if (ia(iCoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                        dble(ia(iCoor)*(ia(iCoor)-1))* &
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
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-rKappa(iZeta)* &
                                                          dble(ia(iCoor)*(ia(iCoor)-1))*Four*Beta(iZeta)*dble(lb)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(iCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*(ia(iCoor)-1)*ib(iCoor)*(ib(iCoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*(ia(iCoor)-1)*ib(jCoor)*(ib(jCoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kCoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*(ia(iCoor)-1)*ib(kCoor)*(ib(kCoor)-1))* &
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
              final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = rKappa(iZeta)* &
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
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                        (Two*Alpha(iZeta))**2*Four*Beta(iZeta)*dble(lb)* &
                                                        Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                        Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(icoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                        (Two*Alpha(iZeta))**2*dble(ib(icoor)*(ib(icoor)-1))* &
                                                        Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                        Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(jcoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                        (Two*Alpha(iZeta))**2*dble(ib(jcoor)*(ib(jcoor)-1))* &
                                                        Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                        Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
              end if
              if (ib(kcoor) > 1) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                        (Two*Alpha(iZeta))**2*dble(ib(kcoor)*(ib(kcoor)-1))* &
                                                        Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                        Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                        Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
              end if
              if (ia(icoor) > 0) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                        dble(ia(icoor))*Two*Alpha(iZeta)* &
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
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(iCoor))*Four*Beta(iZeta)*dble(lb)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(iCoor)*ib(icoor)*(ib(icoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(iCoor)*ib(jcoor)*(ib(jcoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(iCoor)*ib(kcoor)*(ib(kcoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if (ia(jcoor) > 0) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                        dble(ia(jcoor))*Two*Alpha(iZeta)* &
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
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(jCoor))*Four*Beta(iZeta)*dble(lb)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(jCoor)*ib(icoor)*(ib(icoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(jCoor)*ib(jcoor)*(ib(jcoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          Two*Alpha(iZeta)*dble(ia(jCoor))*dble(ib(kcoor)*(ib(kcoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2)
                end if
              end if
              if ((ia(iCoor) > 0) .and. (ia(jCoor) > 0)) then
                final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                        dble(ia(iCoor)*ia(jCoor))* &
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
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-rKappa(iZeta)* &
                                                          dble(ia(iCoor)*ia(jCoor))*Four*Beta(iZeta)*dble(lb)* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(icoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*ia(jCoor)*ib(icoor)*(ib(icoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(jcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*ia(jCoor)*ib(jcoor)*(ib(jcoor)-1))* &
                                                          Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                          Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                          Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
                end if
                if (ib(kcoor) > 1) then
                  final(iZeta,ipa,ipb,I(6,iCoor,jCoor)) = final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+rKappa(iZeta)* &
                                                          dble(ia(iCoor)*ia(jCoor)*ib(kcoor)*(ib(kcoor)-1))* &
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
!  call RecPrt(' S(1)',' ',Final,nDAO,21)
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

            ps = dble(iChTbl(iIrrep,nOp(2))**(iCnt+jCnt))

            ! Get the transf. character of the diff. operator

            ps = ps*dble(iPrmt(nOp(2),iCh))

            ! If the over triangular diff. are needed multiply by two instead

            if ((iCnt /= jCnt) .and. (iCar == jCar) .and. (abs(indgrd(iCar-1,iCnt,iIrrep)) == abs(indgrd(jCar-1,jCnt,iIrrep)))) then
              ps = ps*Two
            end if
            iHess = abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
            Fact = dble(iStab(iCnt)*iStab(jCnt))/dble(nIrrep**2)
            Fact = Fact*ps
            if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) > 0) then
              oj = DDot_(nDAO,DAO,1,final(1,1,1,I(6,iCar,jCar)),1)
              Hess(iHess) = Hess(iHess)+Fact*oj
            else
              Fact = Fact*dble((-1)**(icnt+jcnt))
              oj = DDot_(nDAO,DAO,1,final(1,1,1,I(6,max(iCar,jCar),min(iCar,jCar))),1)
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
