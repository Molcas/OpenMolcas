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
! Copyright (C) 1990,2020, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,ABeq)
!***********************************************************************
!                                                                      *
! Object: to compile the value of the angular part of a basis function *
!         evaluated at a root of the quadrature.                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: One

implicit none
integer, intent(In) :: nZeta, nHer, na
integer iHer, iCar, ia
real*8, intent(In) :: Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer)
real*8, intent(Out) :: Axyz(nZeta,3,nHer,0:na)
#ifdef _DEBUGPRINT_
character*80 Label
#endif
logical ABeq(3)

if (na < 0) then
  call WarningMessage(2,'CrtCmp: na < 0')
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(Label,'(A)') ' In CrtCmp: Axyz(in)'
call RecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1))
call RecPrt(' In CrtCmp: HerR',' ',HerR,1,nHer)
call RecPrt(' In CrtCmp: Zeta',' ',Zeta,nZeta,1)
call RecPrt(' In CrtCmp: A   ',' ',A,1,3)
call RecPrt(' In CrtCmp: P   ',' ',P,nZeta,3)
#endif
Axyz(:,:,:,0) = One
if (na /= 0) then

  do iHer=1,nHer
    do iCar=1,3

      if (ABeq(iCar)) then
        Axyz(:,iCar,iHer,1) = HerR(iHer)*1/sqrt(Zeta(:))
      else
        Axyz(:,iCar,iHer,1) = HerR(iHer)*1/sqrt(Zeta(:))+P(:,iCar)-A(iCar)
      end if

      do ia=2,na
        Axyz(:,iCar,iHer,ia) = Axyz(:,iCar,iHer,1)*Axyz(:,iCar,iHer,ia-1)
      end do

    end do
  end do

end if

#ifdef _DEBUGPRINT_
write(Label,'(A)') ' In CrtCmp: Axyz(out) '
call RecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1))
#endif

return

end subroutine CrtCmp
