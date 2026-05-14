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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine vAssmbl(Rnxyz,Axyz,la,Rxyz,lr,Bxyz,lb,nZeta,HerW,nHer,Temp)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the multipole moment *
!         matrix within the framework of the Gauss-Hermite quadrature. *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: la, lr, lb, nZeta, nHer
real(kind=wp), intent(out) :: Rnxyz(nZeta*3,0:la,0:lb,0:lr), Temp(nZeta*3,nHer)
real(kind=wp), intent(in) :: Axyz(nZeta*3,nHer,0:la), Rxyz(nZeta*3,nHer,0:lr), Bxyz(nZeta*3,nHer,0:lb), HerW(nHer)
integer(kind=iwp) :: ia, ib, iHer, ir
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif

#ifdef _DEBUGPRINT_
call RecPrt(' In vAssmbl:HerW',' ',HerW,1,nHer)
call RecPrt(' In vAssmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1))
call RecPrt(' In vAssmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1))
call RecPrt(' In vAssmbl:Rxyz',' ',Rxyz,nZeta*3,nHer*(lr+1))
#endif

Rnxyz(:,:,:,:) = Zero
do ia=0,la
  do ib=0,lb
    do iHer=1,nHer
      Temp(:,iHer) = Axyz(:,iHer,ia)*Bxyz(:,iHer,ib)*HerW(iHer)
    end do
    do ir=0,lr

      ! Generate the cartesian components of the multipole moment
      ! matrix as a sum of the value of the integrand, evaluated
      ! at a root, times a weight.

      do iHer=1,nHer
        Rnxyz(:,ia,ib,ir) = Rnxyz(:,ia,ib,ir)+Temp(:,iHer)*Rxyz(:,iHer,ir)
      end do

#     ifdef _DEBUGPRINT_
      write(Label,'(A,I2,A,I2,A,I2,A)') ' In vAssmbl: Rnxyz(',ia,',',ib,',',ir,')'
      call RecPrt(Label,' ',Rnxyz(1,ia,ib,ir),nZeta,3)
#     endif
    end do
  end do
end do

end subroutine vAssmbl
