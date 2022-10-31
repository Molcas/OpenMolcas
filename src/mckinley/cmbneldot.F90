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
! Copyright (C) 1991,1992,1995, Roland Lindh                           *
!               1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnEldot(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,Fact,Temp,Alpha,Beta,DAO,iStb,jStb,nOp,rOut,indgrd)
!***********************************************************************
!                                                                      *
! Object: to compute gradient integrals for SC Reaction Fields         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!             Modified for gradient calculations May '95               *
!             Modified for trans. prob. calculations Oct '97           *
!             by Anders Bernhardsson                                   *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri3_Elem, nTri_Elem1
use Symmetry_Info, only: iChBas, iChTbl, nIrrep
use Constants, only: Two, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp, iStb, jStb, nOp(2), indgrd(2,3,3,0:7)
real(kind=wp), intent(in) :: Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), Zeta(nZeta), rKappa(nZeta), Alpha(nZeta), Beta(nZeta), &
                             DAO(nZeta,nTri_Elem1(la),nTri_Elem1(lb))
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp,6), Fact(nZeta), Temp(nZeta)
real(kind=wp), intent(inout) :: rOut(*)
integer(kind=iwp) :: i1, iCar, iCnt, iComp, ihess, iIrrep, ipa, ipb, ir, ix, ixa, ixb, iy, iya, iyaMax, iyb, iybMax, iz, iza, izb, &
                     jCar, nDAO
real(kind=wp) :: Fct, ps, rtemp, xa, xb, ya, yb, za, zb
integer(kind=iwp), external :: iPrmt
real(kind=wp), external :: DDot_

Fact(:) = rKappa*Zeta**(-OneHalf)

! Loop over angular components of the basis set

do ixa=0,la
  iyaMax = la-ixa
  do ixb=0,lb
    iybMax = lb-ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)

        ! Combine multipole moment integrals

        do ix=0,lr
          do iy=0,lr-ix
            if (ixa > 0) then
              xa = -real(ixa,kind=wp)
              Temp(:) = Fact(:)*(Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb,ix)+xa*Rnxyz(:,1,ixa-1,ixb,ix))*Rnxyz(:,2,iya,iyb,iy)
            else
              Temp(:) = Fact(:)*Two*Alpha(:)*Rnxyz(:,1,ixa+1,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              rFinal(:,ipa,ipb,iComp,1) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (ixb > 0) then
              xb = -real(ixb,kind=wp)
              Temp(:) = Fact(:)*(Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1,ix)+xb*Rnxyz(:,1,ixa,ixb-1,ix))*Rnxyz(:,2,iya,iyb,iy)
            else
              Temp(:) = Fact(:)*Two*Beta(:)*Rnxyz(:,1,ixa,ixb+1,ix)*Rnxyz(:,2,iya,iyb,iy)
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              rFinal(:,ipa,ipb,iComp,4) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (iya > 0) then
              ya = -real(iya,kind=wp)
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*(Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb,iy)+ya*Rnxyz(:,2,iya-1,iyb,iy))
            else
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Two*Alpha(:)*Rnxyz(:,2,iya+1,iyb,iy)
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              rFinal(:,ipa,ipb,iComp,2) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (iyb > 0) then
              yb = -real(iyb,kind=wp)
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*(Two*Beta(:)*Rnxyz(:,2,iya,iyb+1,iy)+yb*Rnxyz(:,2,iya,iyb-1,iy))
            else
              Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Two*Beta(:)*Rnxyz(:,2,iya,iyb+1,iy)
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              rFinal(:,ipa,ipb,iComp,5) = Temp(:)*Rnxyz(:,3,iza,izb,iz)
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              if (iza > 0) then
                za = -real(iza,kind=wp)
                rFinal(:,ipa,ipb,iComp,3) = Temp(:)*(Two*Alpha(:)*Rnxyz(:,3,iza+1,izb,iz)+za*Rnxyz(:,3,iza-1,izb,iz))
              else
                rFinal(:,ipa,ipb,iComp,3) = Temp(:)*Two*Alpha(:)*Rnxyz(:,3,iza+1,izb,iz)
              end if
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            Temp(:) = Fact(:)*Rnxyz(:,1,ixa,ixb,ix)*Rnxyz(:,2,iya,iyb,iy)

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = C_Ind(ir,ix,iz)+nTri3_Elem(ir)
              if (izb > 0) then
                zb = -real(izb,kind=wp)
                rFinal(:,ipa,ipb,iComp,6) = Temp(:)*(Two*Beta(:)*Rnxyz(:,3,iza,izb+1,iz)+zb*Rnxyz(:,3,iza,izb-1,iz))
              else
                rFinal(:,ipa,ipb,iComp,6) = Temp(:)*Two*Beta(:)*Rnxyz(:,3,iza,izb+1,iz)
              end if
            end do
          end do
        end do

      end do
    end do
  end do
end do

nDAO = nZeta*nTri_Elem1(la)*nTri_Elem1(lb)
do iIrrep=0,nIrrep-1
  do iCnt=1,2
    do iCar=1,3
      do jCar=1,3
        icomp = jcar+1
        if (iCnt == 1) then
          i1 = iCar
          ps = real(iChTbl(iIrrep,nOp(1)),kind=wp)
          ps = ps*real(iPrmt(nOp(1),iChBas(1+iCar)),kind=wp)
          Fct = real(iStb,kind=wp)/real(nIrrep,kind=wp)
        else
          i1 = iCar+3
          ps = real(iChTbl(iIrrep,nOp(2)),kind=wp)
          ps = ps*real(iPrmt(nOp(2),iChBas(1+iCar)),kind=wp)
          Fct = ps*real(jStb,kind=wp)/real(nIrrep,kind=wp)
        end if

        if (IndGrd(iCnt,iCar,jCar,iIrrep) /= 0) then
          ihess = indgrd(icnt,icar,jcar,iirrep)
          rtemp = DDot_(nDAO,DAO,1,rFinal(:,:,:,icomp,i1),1)
          rOut(iHess) = rOut(iHess)+Fct*rtemp
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnEldot
