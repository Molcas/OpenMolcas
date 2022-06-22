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
!***********************************************************************

subroutine Contact(Zeta,P,nZeta,A,Axyz,la,RB,Bxyz,lb,Ccoor,lOper,iChO,nIC,Array,final,iStabM,nStabM,nComp,rKappa)
!***********************************************************************
!                                                                      *
! Object: to compoute the 1-electron contact term.                     *
!                                                                      *
!     Author: Roland Lindh, Dept. Of Theoretical Chemistry,            *
!             University of Lund, Sweden, February '91                 *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC), rKappa(nZeta), Ccoor(3), Array(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), &
       Axyz(nZeta,3,0:la), Bxyz(nZeta,3,0:lb), Zeta(nZeta), P(nZeta,3), A(3), RB(3), TC(3)
integer iStabM(0:nStabM-1), iStabO(0:7), iDCRT(0:7), lOper(nComp), iChO(nComp)
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 170
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In Contact: rKappa',' ',rKappa,nZeta,1)
  call RecPrt(' In Contact: Zeta',' ',Zeta,nZeta,1)
  call RecPrt(' In Contact: P',' ',P,nZeta,3)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
llOper = lOper(1)
do iComp=2,nComp
  llOper = ior(llOper,lOper(iComp))
end do
call SOS(iStabO,nStabO,llOper)
call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

do lDCRT=0,nDCRT-1
  call OA(iDCRT(lDCRT),Ccoor,TC)

  call dcopy_(nZeta*nElem(la)*nElem(lb),[Zero],0,Array,1)

  ! Compute the value of the angular components associated
  ! to the basis functions centered on the first center.

  call dcopy_(nZeta*3,[One],0,Axyz(1,1,0),1)

  if (la /= 0) then
    do iCar=1,3

      do iZeta=1,nZeta
        Axyz(iZeta,iCar,1) = TC(iCar)-A(iCar)
      end do

      do ia=2,la
        do iZeta=1,nZeta
          Axyz(iZeta,iCar,ia) = Axyz(iZeta,iCar,1)*Axyz(iZeta,iCar,ia-1)
        end do
      end do

    end do
  end if

  ! Compute the value of the angular components associated to
  ! the basis functions centered on the second center.

  call dcopy_(nZeta*3,[One],0,Bxyz(1,1,0),1)

  ! Modify z-component to carry the the exponential contribution.

  do iZeta=1,nZeta
    Bxyz(iZeta,3,0) = exp(-Zeta(iZeta)*((TC(1)-P(iZeta,1))**2+(TC(2)-P(iZeta,2))**2+(TC(3)-P(iZeta,3))**2))
  end do

  if (lb /= 0) then
    do iCar=1,3

      do iZeta=1,nZeta
        Bxyz(iZeta,iCar,1) = TC(iCar)-RB(iCar)
      end do

      do ib=2,lb
        do iZeta=1,nZeta
          Bxyz(iZeta,iCar,ib) = Bxyz(iZeta,iCar,1)*Bxyz(iZeta,iCar,ib-1)
        end do
      end do
    end do

    ! Modify z-components with the exponential contribution

    do ib=1,lb
      do iZeta=1,nZeta
        Bxyz(iZeta,3,ib) = Bxyz(iZeta,3,ib)*Bxyz(iZeta,3,0)
      end do
    end do
  end if

  ! Combine contributions from the various angular components.

  do ixa=la,0,-1
    do ixb=lb,0,-1
      do iya=la-ixa,0,-1
        iza = la-ixa-iya
        ipa = Ind(la,ixa,iza)
        do iyb=lb-ixb,0,-1
          izb = lb-ixb-iyb
          ipb = Ind(lb,ixb,izb)
          do iZeta=1,nZeta
            Array(iZeta,ipa,ipb) = Array(iZeta,ipa,ipb)+rKappa(iZeta)*Axyz(iZeta,1,ixa)*Axyz(iZeta,2,iya)*Axyz(iZeta,3,iza)* &
                                   Bxyz(iZeta,1,ixb)*Bxyz(iZeta,2,iyb)*Bxyz(iZeta,3,izb)
          end do
        end do
      end do
    end do
  end do

  ! Accumulate contributions

  nOp = NrOpr(iDCRT(lDCRT))
  call SymAdO(Array,nZeta,la,lb,nComp,final,nIC,nOp,lOper,iChO,One)

end do

return

end subroutine Contact
