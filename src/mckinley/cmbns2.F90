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

subroutine CmbnS2(Rnxyz,nZeta,la,lb,Zeta,rKappa,rFinal,Alpha,Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
!***********************************************************************
!                                                                      *
! Object: compute the 2nd derivative of the overlap matrix.            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: C_Ind, iTri, nTri_Elem1
use Symmetry_Info, only: iChTbl, nIrrep
use Constants, only: One, Two, Four, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nHess, IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1), indgrd(0:2,0:1,0:nirrep-1), iu, iv, &
                                 nOp(2)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+2,0:lb+2), Zeta(nZeta), Alpha(nZeta), DAO(nZeta,nTri_Elem1(la),nTri_Elem1(lb))
real(kind=wp), intent(inout) :: rKappa(nZeta), Hess(nHess)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6)
logical(kind=iwp), intent(in) :: IfHss(0:1,0:2,0:1,0:2)
integer(kind=iwp) :: i, ia(3), iax, iay, ib(3), ibx, iby, iCar, iCh, iCnt, iCoor, iHess, iIrrep, ipa, ipb, istb(0:1), iStop, &
                     iyaMax, iybMax, jCar, jCnt, jCoor, kCoor, nDAO
real(kind=wp) :: Fact, ps, rtemp
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_

!EAW 970912 ixyz = ixLoc(DAO(1,1,1))
!iRout = 134
!iPrint = nPrint(iRout)
iStb(0) = iu
iStb(1) = iv

rKappa(:) = rKappa*Zeta**(-OneHalf)
!if (iPrint >= 99) then
!  call RecPrt(' In CmbnS2: Zeta  ',' ',Zeta,1,nZeta)
!  call RecPrt(' In CmbnS2: rKappa',' ',rKappa,1,nZeta)
!  call RecPrt(' In CmbnS2: Alpha ',' ',Alpha,1,nZeta)
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
          if (IfHss(0,iCoor-1,0,iCoor-1)) then
            rFinal(:,ipa,ipb,i) = rKappa(:)*((Two*Alpha(:))**2*Rnxyz(:,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                               Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                               Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))- &
                                             Two*Alpha(:)*Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                          Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                          Rnxyz(:,kCoor,ia(kCoor),ib(kCoor)))
            if (ia(iCoor) > 0) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)-rKappa(:)*Four*Alpha(:)*real(ia(iCoor),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor),ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
            if (ia(iCoor) > 1) rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)+rKappa(:)*real(ia(iCoor)*(ia(iCoor)-1),kind=wp)* &
                                                     Rnxyz(:,iCoor,ia(iCoor)-2,ib(iCoor))* &
                                                     Rnxyz(:,jCoor,ia(jCoor),ib(jCoor))* &
                                                     Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
          end if
        end do

        ! Integrals like dI/dxdz

        do iCoor=2,3
          do jCoor=1,iCoor-1
            i = iTri(iCoor,jCoor)
            if (IfHss(0,iCoor-1,0,jCoor-1)) then
              rFinal(:,ipa,ipb,i) = rKappa
              do kCoor=1,3
                if ((kCoor == iCoor) .or. (kCoor == jCoor)) then
                  if (ia(kCoor) > 0) then
                    rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)*(Two*Alpha(:)*Rnxyz(:,kCoor,ia(kCoor)+1,ib(kCoor))- &
                                                               real(ia(kCoor),kind=wp)*Rnxyz(:,kCoor,ia(kCoor)-1,ib(kCoor)))
                  else
                    rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)*Two*Alpha(:)*Rnxyz(:,kCoor,ia(kCoor)+1,ib(kCoor))
                  end if
                else
                  rFinal(:,ipa,ipb,i) = rFinal(:,ipa,ipb,i)*Rnxyz(:,kCoor,ia(kCoor),ib(kCoor))
                end if
              end do
            end if
          end do
        end do
      end do
    end do
  end do
end do

! Trace the Hessian integrals
nDAO = nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
!if (iPrint >= 99) then
!  call RecPrt(' S(1)',' ',rFinal,nDAO,21)
!  call RecPrt('   D ','(6f12.6)',DAO(1,1,1),nDAO,1)
!end if
do iIrrep=0,nIrrep-1
  do iCnt=0,1
    do iCar=1,3
      do jCnt=0,1
        if (iCnt == jCnt) then
          iStop = iCar
        else
          iStop = 3
        end if
        do jCar=1,iStop
          i = iTri(iCar,jCar)
          if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) /= 0) then

            ! Accumulate contribution to the Hessian

            ! Get the characteristics of the diff operator

            iCh = ieor(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)

            ! Get the character of the operator in the present irrep

            ps = real(iChTbl(iIrrep,nOp(2))**(iCnt+jCnt),kind=wp)

            ! Get the transf. character of the diff. operator

            ps = ps*real(iPrmt(nOp(2),iCh),kind=wp)

            ! If the over triangular diff. are needed  multiply by two instead
            ! Because of that x2x1 y2y1 z2z1 just appear ones in the (1,2)
            ! "subhessian".

            if ((iCnt /= jCnt) .and. (iCar == jCar) .and. &
                (abs(indgrd(iCar-1,iCnt,iIrrep)) == abs(indgrd(jCar-1,jCnt,iIrrep)))) ps = ps*Two
            iHess = abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
            Fact = real(iStb(iCnt)*iStb(jCnt),kind=wp)/real(nIrrep**2,kind=wp)
            Fact = Fact*ps
            if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) > 0) then
              rtemp = DDot_(nDAO,DAO,1,rFinal(:,:,:,i),1)
              Hess(iHess) = Hess(iHess)+Fact*rtemp
            else
              Fact = Fact*(-One)**(iCnt+jCnt)
              rtemp = DDot_(nDAO,DAO,1,rFinal(:,:,:,i),1)
              Hess(iHess) = Hess(iHess)+Fact*rtemp
            end if
          end if
        end do
      end do
    end do
  end do
end do

return

end subroutine CmbnS2
