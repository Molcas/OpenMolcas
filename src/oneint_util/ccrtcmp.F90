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
! Copyright (C) 1990,2015, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine CCrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,kVector)
!***********************************************************************
!                                                                      *
! Object: to compile the value of the angular part of a basis function *
!         evaluated at a root of the quadrature.                       *
!                                                                      *
! Called from: PrpInt                                                  *
!                                                                      *
! Calling    :                                                         *
!              RecPrt                                                  *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!                                                                      *
!             Roland Lindh, Uppsala Universitet, Uppsala Sweden        *
!             December 2015.                                           *
!             Modification to wave vectors and complex representation. *
!***********************************************************************

use Constants, only: Two, cOne
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, na, nHer
real(kind=wp), intent(in) :: Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer), kVector(3)
complex(kind=wp), intent(out) :: Axyz(nZeta,3,nHer,0:na)
integer(kind=iwp) :: ia, iCar, iHer
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif

if (na < 0) then
  call WarningMessage(2,'CCrtCmp: na < 0')
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In CCrtCmp: HerR',' ',HerR,1,nHer)
call RecPrt(' In CCrtCmp: Zeta',' ',Zeta,nZeta,1)
call RecPrt(' In CCrtCmp: A   ',' ',A,1,3)
call RecPrt(' In CCrtCmp: P   ',' ',P,nZeta,3)
call RecPrt(' In CCrtCmp: KVec',' ',kVector,1,3)
#endif

Axyz(:,:,:,0) = cOne

if (na /= 0) then
  do iHer=1,nHer
    do iCar=1,3

      Axyz(:,iCar,iHer,1) = cmplx(HerR(iHer)/sqrt(Zeta)+P(:,iCar)-A(iCar),kVector(iCar)/(Two*Zeta),kind=wp)

      do ia=2,na
        Axyz(:,iCar,iHer,ia) = Axyz(:,iCar,iHer,1)*Axyz(:,iCar,iHer,ia-1)
      end do

    end do
  end do
end if

#ifdef _DEBUGPRINT_
write(Label,'(A)') ' In CCrtCmp: Axyz '
call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'R')
call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'I')
#endif

end subroutine CCrtCmp
