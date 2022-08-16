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

subroutine CmbnS2(Rnxyz,nZeta,la,lb,Zeta,rKappa,final,Alpha,Beta,Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
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
integer IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1), istb(0:1), nOp(2), ia(3), ib(3), indgrd(0:2,0:1,0:nirrep-1)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
! Index in the triang. local hessian
I(i1,i2) = i1*(i1-1)/2+i2

!EAW 970912 ixyz = ixLoc(DAO(1,1,1))
!iRout = 134
!iPrint = nPrint(iRout)
iStb(0) = iu
iStb(1) = iv
!call GetMem(' Enter CmbnS2','LIST','REAL',iDum,iDum)

exp32 = -Three/Two
do iZeta=1,nZeta
  rKappa(iZeta) = rKappa(iZeta)*Zeta(iZeta)**exp32
end do
!if (iPrint >= 99) then
!  call RecPrt(' In CmbnS2: Zeta  ',' ',Zeta,1,nZeta)
!  call RecPrt(' In CmbnS2: rKappa',' ',rKappa,1,nZeta)
!  call RecPrt(' In CmbnS2: Alpha ',' ',Alpha,1,nZeta)
!  call RecPrt(' In CmbnS2: Beta  ',' ',Beta,1,nZeta)
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
              final(iZeta,ipa,ipb,I(iCoor,iCoor)) = rKappa(iZeta)*((Two*Alpha(iZeta))**2*Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))* &
                                                                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))- &
                                                                   Two*Alpha(iZeta)*Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              if (ia(iCoor) > 0) then
                final(iZeta,ipa,ipb,I(iCoor,iCoor)) = final(iZeta,ipa,ipb,I(iCoor,iCoor))- &
                                                      rKappa(iZeta)*(Four*Alpha(iZeta)*dble(ia(iCoor))* &
                                                                     Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))* &
                                                                     Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))* &
                                                                     Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
              end if
              if (ia(iCoor) > 1) then
                final(iZeta,ipa,ipb,I(iCoor,iCoor)) = final(iZeta,ipa,ipb,I(iCoor,iCoor))+ &
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
            if (IfHss(0,iCoor-1,0,jCoor-1)) then
              do kCoor=1,3
                do iZeta=1,nZeta
                  if (kCoor == 1) then
                    final(iZeta,ipa,ipb,I(iCoor,jCoor)) = rKappa(iZeta)
                  end if
                  if ((kCoor == iCoor) .or. (kCoor == jCoor)) then
                    rIc = Two*Alpha(iZeta)*Rnxyz(iZeta,kCoor,ia(kCoor)+1,ib(kCoor))

                    if (ia(kCoor) > 0) rIc = rIc-dble(ia(kCoor))*Rnxyz(iZeta,kCoor,ia(kCoor)-1,ib(kCoor))

                    final(iZeta,ipa,ipb,I(iCoor,jCoor)) = final(iZeta,ipa,ipb,I(iCoor,jCoor))*rIc
                  else
                    final(iZeta,ipa,ipb,I(iCoor,jCoor)) = final(iZeta,ipa,ipb,I(iCoor,jCoor))*Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor))
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

! Trace the Hessian integrals
nDAO = nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
!if (iPrint >= 99) then
!  call RecPrt(' S(1)',' ',Final,nDAO,21)
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
          if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) /= 0) then

            ! Accumulate contribution to the Hessian

            ! Get the characteristics of the diff operator

            iCh = ieor(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)

            ! Get the character of the operator in the present irrep

            ps = dble(iChTbl(iIrrep,nOp(2))**(iCnt+jCnt))

            ! Get the transf. character of the diff. operator

            ps = ps*dble(iPrmt(nOp(2),iCh))

            ! If the over triangular diff. are needed  multiply by two instead
            ! Because of that x2x1 y2y1 z2z1 just appear ones in the (1,2)
            ! "subhessian".

            if ((iCnt /= jCnt) .and. (iCar == jCar) .and. (abs(indgrd(iCar-1,iCnt,iIrrep)) == abs(indgrd(jCar-1,jCnt,iIrrep)))) then
              ps = ps*Two
            end if
            iHess = abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
            Fact = dble(iStb(iCnt)*iStb(jCnt))/dble(nIrrep**2)
            Fact = Fact*ps
            if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep) > 0) then
              rtemp = DDot_(nDAO,DAO,1,final(1,1,1,I(max(iCar,jCar),min(iCar,jCar))),1)
              Hess(iHess) = Hess(iHess)+Fact*rtemp
            else
              Fact = Fact*dble((-1)**(icnt+jcnt))
              rtemp = DDot_(nDAO,DAO,1,final(1,1,1,I(max(iCar,jCar),min(iCar,jCar))),1)
              Hess(iHess) = Hess(iHess)+Fact*rtemp
            end if
          end if
        end do
      end do
    end do
  end do
end do

!call GetMem(' Exit CmbnS2','LIST','REAL',iDum,iDum)

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(Beta)

end subroutine CmbnS2
