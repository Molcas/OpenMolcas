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
subroutine vCrtCmp(Zeta12,P,nZeta,A,Axyz,na,HerR,nHer,ABeq)
!***********************************************************************
!                                                                      *
! Object: to compile the value of the angular part of a basis function *
!         evaluated at a root of the quadrature.                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, na, nHer
real(kind=wp), intent(in) :: Zeta12(nZeta), P(nZeta,3), A(3), HerR(nHer)
real(kind=wp), intent(out) :: Axyz(nZeta,3,nHer,0:na)
logical(kind=iwp), intent(in) :: ABeq(3)
integer(kind=iwp) :: ia, iCar, iHer
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif

#ifdef _DEBUGPRINT_
call RecPrt(' In vCrtCmp: HerR',' ',HerR,1,nHer)
call RecPrt(' In vCrtCmp: Zeta',' ',Zeta12,nZeta,1)
call RecPrt(' In vCrtCmp: A   ',' ',A,1,3)
call RecPrt(' In vCrtCmp: P   ',' ',P,nZeta,3)
#endif
Axyz(:,:,:,0) = One
if (na /= 0) then

  do iHer=1,nHer
    do iCar=1,3

      if (ABeq(iCar)) then
        Axyz(:,iCar,iHer,1) = HerR(iHer)*Zeta12(:)
      else
        Axyz(:,iCar,iHer,1) = HerR(iHer)*Zeta12(:)+P(:,iCar)-A(iCar)
      end if

      do ia=2,na
        Axyz(:,iCar,iHer,ia) = Axyz(:,iCar,iHer,1)*Axyz(:,iCar,iHer,ia-1)
      end do

    end do
  end do

end if

#ifdef _DEBUGPRINT_
do ia=0,na
  do iHer=1,nHer
    write(Label,'(A,I2,A,I2,A)') ' In vCrtCmp: Axyz (iHer=',iHer,',ia=',ia,')'
    call RecPrt(Label,' ',Axyz(:,:,iHer,ia),nZeta,3)
  end do
end do
#endif

end subroutine vCrtCmp
