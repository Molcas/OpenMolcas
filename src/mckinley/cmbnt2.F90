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

use Index_Functions, only: C_Ind, iTri, nTri_Elem1
use Symmetry_Info, only: nIrrep, iChTbl
use Constants, only: One, Two, Four, Six, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nHess, IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1), indgrd(0:2,0:1,0:nirrep-1), iu, iv, &
                                 nOp(2)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+2,0:lb+2), Zeta(nZeta), Alpha(nZeta), Beta(nZeta), &
                             DAO(nZeta,nTri_Elem1(la),nTri_Elem1(lb))
real(kind=wp), intent(inout) :: rKappa(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6), Hess(nHess)
logical(kind=iwp), intent(in) :: IfHss(0:1,0:2,0:1,0:2)
integer(kind=iwp) :: i, ia(3), iax, iay, ib(3), ibx, iby, iCar, ich, iCnt, iCoor, iHess, iIrrep, iMax, ipa, ipb, istab(0:1), &
                     iyaMax, iybMax, jCar, jCnt, jCoor, kCoor, nDAO
real(kind=wp) :: Fact, oj, ps
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_

!iRout = 134
iStab(0) = iu
iStab(1) = iv
!iPrint = nPrint(iRout)

rKappa(:) = Half*rKappa*Zeta**(-OneHalf)
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
          i = iTri(iCoor,iCoor)
          if (.not. IfHss(0,iCoor-1,0,iCoor-1)) cycle
          rFinal(:,ipa,ipb,i) = rKappa(:)*((Two*Alpha(:))**2*((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor)+2)* &
                                                                                Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                                Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                                Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                                                Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                                Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                                Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                              Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                          Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                          Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))- &
                                           Two*Alpha(:)*((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor),ib(iCoor)+2)* &
                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                           Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                           Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                         Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))))
          if (lb > 0) &
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                  ((Two*Alpha(:))**2*Four*Beta(:)*real(lb,kind=wp)*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                                   Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                   Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))- &
                                   Two*Alpha(:)*(Four*Beta(:))*real(lb,kind=wp)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                                Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
          if (ib(icoor) > 1) &
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                  ((Two*Alpha(:))**2*real(ib(icoor)*(ib(icoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor)-2)* &
                                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))- &
                                   Two*Alpha(:)*real(ib(icoor)*(ib(icoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                                                      Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                      Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
          if (ib(jcoor) > 1) &
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                  ((Two*Alpha(:))**2*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))- &
                                   Two*Alpha(:)*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                                      Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                                                      Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
          if (ib(kcoor) > 1) &
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                  ((Two*Alpha(:))**2*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                                           Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                           Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)- &
                                   Two*Alpha(:)*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                                                      Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                                                      Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2))
          if (ia(iCoor) > 0) then
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                  ((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor),ib(iCoor)+2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                   Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                               Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                              Four*Beta(:)*real(lb,kind=wp)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                              Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                              Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                              Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(iCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     real(ib(icoor)*(ib(icoor)-1),kind=wp)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor)-2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(jCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     real(ib(jcoor)*(ib(jcoor)-1),kind=wp)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(kCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     real(ib(kcoor)*(ib(kcoor)-1),kind=wp)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          end if
          if (ia(iCoor) > 1) then
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)*real(ia(iCoor)*(ia(iCoor)-1),kind=wp)* &
                                  ((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor)+2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)+2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                   Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                               Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                              real(ia(iCoor)*(ia(iCoor)-1),kind=wp)*Four*Beta(:)*real(lb,kind=wp)* &
                                              Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                              Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                              Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(iCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*(ia(iCoor)-1)*ib(iCoor)*(ib(iCoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor)-2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(jCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*(ia(iCoor)-1)*ib(jCoor)*(ib(jCoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor)-2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(kCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*(ia(iCoor)-1)*ib(kCoor)*(ib(kCoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          end if
        end do

        ! Integrals like dI/dxdz

        do kCoor=1,3
          iCoor = mod(kCoor,3)+1
          jCoor = mod(iCoor,3)+1
          iMax = max(iCoor,jCoor)
          jCoor = min(iCoor,jCoor)
          iCoor = iMax
          i = iTri(iCoor,jCoor)
          if (.not. IfHss(0,iCoor-1,0,jCoor-1)) cycle
          rFinal(:,ipa,ipb,i) = rKappa(:)*(Two*Alpha(:))**2*((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor)+2)* &
                                                                               Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                               Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                                               Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor)+2)* &
                                                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                                               Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                                               Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                                             Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                                         Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                                         Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
          if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)*(Two*Alpha(:))**2*Four*Beta(:)*real(lb,kind=wp)* &
                                            Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                            Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                            Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
          if (ib(icoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                   (Two*Alpha(:))**2*real(ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                   Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                   Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                   Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
          if (ib(jcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                   (Two*Alpha(:))**2*real(ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                   Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                   Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                   Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
          if (ib(kcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                   (Two*Alpha(:))**2*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                   Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                   Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                   Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          if (ia(icoor) > 0) then
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)*real(ia(icoor),kind=wp)*Two*Alpha(:)* &
                                  ((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor)+2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor)+2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                   Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                               Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                              Two*Alpha(:)*real(ia(iCoor))*Four*Beta(:)*real(lb,kind=wp)* &
                                              Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                              Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                              Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(icoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(iCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(jcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(iCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor)-2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(kcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(iCoor)*ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)+1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          end if
          if (ia(jcoor) > 0) then
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)*real(ia(jcoor),kind=wp)*Two*Alpha(:)* &
                                  ((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor)+2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor)+2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                   Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                               Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                              Two*Alpha(:)*real(ia(jCoor))*Four*Beta(:)*real(lb,kind=wp)* &
                                              Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                              Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                              Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(icoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(jCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor)-2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(jcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(jCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(kcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                                     Two*Alpha(:)*real(ia(jCoor),kind=wp)*real(ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)+1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          end if
          if ((ia(iCoor) > 0) .and. (ia(jCoor) > 0)) then
            rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)*real(ia(iCoor)*ia(jCoor),kind=wp)* &
                                  ((Two*Beta(:))**2*(Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor)+2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor)+2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))+ &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)+2))- &
                                   Six*Beta(:)*Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                               Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (lb > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)* &
                                              real(ia(iCoor)*ia(jCoor),kind=wp)*Four*Beta(:)*real(lb,kind=wp)* &
                                              Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                              Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                              Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(icoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*ia(jCoor)*ib(icoor)*(ib(icoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor)-2)* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(jcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*ia(jCoor)*ib(jcoor)*(ib(jcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor)-2)* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ib(kcoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)* &
                                                     real(ia(iCoor)*ia(jCoor)*ib(kcoor)*(ib(kcoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-1,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor)-1,ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)-2)
          end if
        end do
      end do
    end do
  end do
end do

! Trace the Hessian integrals

nDAO = nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
!if (iPrint >= 99) then
!  call RecPrt(' S(1)',' ',rFinal,nDAO,21)
!  call RecPrt('   D ',' ',DAO,nDAO,1)
!end if
do iIrrep=0,nIrrep-1
  do iCnt=0,1
    do iCar=1,3
      do jCnt=0,1
        do jCar=1,3
          i = iTri(iCar,jCar)
          if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) /= 0) then

            ! Accumulate contribution to the Hessian

            ! Get the characteristics of the diff operator

            iCh = ieor(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)

            ! Get the character of the operator in the present irrep

            ps = real(iChTbl(iIrrep,nOp(2))**(iCnt+jCnt),kind=wp)

            ! Get the transf. character of the diff. operator

            ps = ps*real(iPrmt(nOp(2),iCh),kind=wp)

            ! If the over triangular diff. are needed multiply by two instead

            if ((iCnt /= jCnt) .and. (iCar == jCar) .and. &
                (abs(indgrd(iCar-1,iCnt,iIrrep)) == abs(indgrd(jCar-1,jCnt,iIrrep)))) ps = ps*Two
            iHess = abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
            Fact = real(iStab(iCnt)*iStab(jCnt),kind=wp)/real(nIrrep**2,kind=wp)
            Fact = Fact*ps
            if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) <= 0) Fact = Fact*(-One)**(iCnt+jCnt)
            oj = DDot_(nDAO,DAO,1,rFinal(:,:,:,i),1)
            Hess(iHess) = Hess(iHess)+Fact*oj
          end if
        end do
      end do
    end do
  end do
end do

return

end subroutine CmbnT2
