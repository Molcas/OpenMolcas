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

subroutine CmbnEldot(Rnxyz,nZeta,la,lb,lr,Zeta,rKappa,final,nComp,Fact,Temp,Alpha,Beta,DAO,iStb,jStb,nOp,rout,indgrd)
!***********************************************************************
!                                                                      *
! Object: to compute gradient integrals for SC Reaction Fields         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Modified for reaction field calculations July '92        *
!             Modified for gradient calculations May '95               *
!             Modified for trans. prob.   calculations Oct '97         *
!             by Anders Bernhardsson                                   *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChTbl, iChBas

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
integer nOp(2), indgrd(2,3,3,0:7)
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp,6), Zeta(nZeta), rKappa(nZeta), Fact(nZeta), Temp(nZeta), &
       Rnxyz(nZeta,3,0:la+1,0:lb+1,0:lr), DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), Alpha(nZeta), Beta(nZeta), rout(*)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6

tTwo = Two

do iZeta=1,nZeta
  Fact(iZeta) = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)
end do

! Loop over angular components of the basis set

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

        ! Combine multipole moment integrals

        do ix=0,lr
          do iy=0,lr-ix
            if (ixa > 0) then
              xa = -ixa
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*(tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)+xa*Rnxyz(iZeta,1,ixa-1,ixb,ix))* &
                              Rnxyz(iZeta,2,iya,iyb,iy)
              end do
            else
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,1,ixa+1,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
              end do
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              do iZeta=1,nZeta
                final(iZeta,ipa,ipb,iComp,1) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
              end do
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (ixb > 0) then
              xb = -ixb
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)+xb*Rnxyz(iZeta,1,ixa,ixb-1,ix))* &
                              Rnxyz(iZeta,2,iya,iyb,iy)
              end do
            else
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*tTwo*Beta(iZeta)*Rnxyz(iZeta,1,ixa,ixb+1,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
              end do
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              do iZeta=1,nZeta
                final(iZeta,ipa,ipb,iComp,4) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
              end do
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (iya > 0) then
              ya = -iya
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*(tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)+ &
                              ya*Rnxyz(iZeta,2,iya-1,iyb,iy))
              end do
            else
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,2,iya+1,iyb,iy)
              end do
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              do iZeta=1,nZeta
                final(iZeta,ipa,ipb,iComp,2) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
              end do
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            if (iyb > 0) then
              yb = -iyb
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)+ &
                              yb*Rnxyz(iZeta,2,iya,iyb-1,iy))
              end do
            else
              do iZeta=1,nZeta
                Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*tTwo*Beta(iZeta)*Rnxyz(iZeta,2,iya,iyb+1,iy)
              end do
            end if

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              do iZeta=1,nZeta
                final(iZeta,ipa,ipb,iComp,5) = Temp(iZeta)*Rnxyz(iZeta,3,iza,izb,iz)
              end do
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            do iZeta=1,nZeta
              Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
            end do

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              if (iza > 0) then
                za = -iza
                do iZeta=1,nZeta
                  final(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)*(tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)+ &
                              za*Rnxyz(iZeta,3,iza-1,izb,iz))
                end do
              else
                do iZeta=1,nZeta
                  final(iZeta,ipa,ipb,iComp,3) = Temp(iZeta)*tTwo*Alpha(iZeta)*Rnxyz(iZeta,3,iza+1,izb,iz)
                end do
              end if
            end do
          end do
        end do
        do ix=0,lr
          do iy=0,lr-ix
            do iZeta=1,nZeta
              Temp(iZeta) = Fact(iZeta)*Rnxyz(iZeta,1,ixa,ixb,ix)*Rnxyz(iZeta,2,iya,iyb,iy)
            end do

            do ir=ix+iy,lr
              iz = ir-ix-iy
              iComp = Ind(ir,ix,iz)+iOff(ir)
              if (izb > 0) then
                zb = -izb
                do iZeta=1,nZeta
                  final(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)*(tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)+ &
                              zb*Rnxyz(iZeta,3,iza,izb-1,iz))
                end do
              else
                do iZeta=1,nZeta
                  final(iZeta,ipa,ipb,iComp,6) = Temp(iZeta)*tTwo*Beta(iZeta)*Rnxyz(iZeta,3,iza,izb+1,iz)
                end do
              end if
            end do
          end do
        end do

      end do
    end do
  end do
end do

nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
do iIrrep=0,nIrrep-1
  do iCnt=1,2
    do iCar=1,3
      do jCar=1,3
        icomp = jcar+1
        if (iCnt == 1) then
          i1 = iCar
          ps = dble(iChTbl(iIrrep,nOp(1)))
          ps = ps*dble(iPrmt(nOp(1),iChBas(1+iCar)))
          Fct = dble(iStb)/dble(nIrrep)
        else
          i1 = iCar+3
          ps = dble(iChTbl(iIrrep,nOp(2)))
          ps = ps*dble(iPrmt(nOp(2),iChBas(1+iCar)))
          Fct = ps*dble(jStb)/dble(nIrrep)
        end if

        if (IndGrd(iCnt,iCar,jCar,iIrrep) /= 0) then
          ihess = indgrd(icnt,icar,jcar,iirrep)
          rtemp = DDot_(nDAO,DAO,1,final(1,1,1,icomp,i1),1)
          rout(iHess) = rOut(iHess)+Fct*rtemp
        end if

      end do
    end do
  end do
end do

return

end subroutine CmbnEldot
