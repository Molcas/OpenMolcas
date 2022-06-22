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

subroutine CCrtCmp(Zeta,P,nZeta,A,Axyz,na,HerR,nHer,KVector)
!***********************************************************************
!                                                                      *
! Object: to compile the value of the angular part of a basis function *
!         evaluated at a root of the quadrature.                       *
!                                                                      *
! Called from: PrpInt                                                  *
!                                                                      *
! Calling    :                                                         *
!              RecPrt                                                  *
!              DCopy (ESSL)                                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!                                                                      *
!             Roland Lindh, Uppsala Universitet, Uppsala Sweden        *
!             December 2015.                                           *
!             Modification to wave vectors and complex representation. *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 Zeta(nZeta), P(nZeta,3), A(3), HerR(nHer), KVector(3)
complex*16 Axyz(nZeta,3,nHer,0:na)
character*80 Label

iRout = 116
iPrint = nPrint(iRout)

if (na < 0) then
  call WarningMessage(2,'CCrtCmp: na < 0')
  call Abend()
end if
if (iPrint >= 99) then
  call RecPrt(' In CCrtCmp: HerR',' ',HerR,1,nHer)
  call RecPrt(' In CCrtCmp: Zeta',' ',Zeta,nZeta,1)
  call RecPrt(' In CCrtCmp: A   ',' ',A,1,3)
  call RecPrt(' In CCrtCmp: P   ',' ',P,nZeta,3)
  call RecPrt(' In CCrtCmp: KVec',' ',KVector,1,3)
end if
do iHer=1,nHer
  do iCar=1,3
    do iZeta=1,nZeta
      Axyz(iZeta,iCar,iHer,0) = DCMPLX(One,Zero)
    end do
  end do
end do

if (na /= 0) then
  do iHer=1,nHer
    do iCar=1,3

      do iZeta=1,nZeta
        Axyz(iZeta,iCar,iHer,1) = DCMPLX(HerR(iHer)*1/sqrt(Zeta(iZeta))+P(iZeta,iCar)-A(iCar),KVector(iCar)/(Two*Zeta(iZeta)))
      end do

      do ia=2,na
        do iZeta=1,nZeta
          Axyz(iZeta,iCar,iHer,ia) = Axyz(iZeta,iCar,iHer,1)*Axyz(iZeta,iCar,iHer,ia-1)
        end do
      end do

    end do
  end do
end if

if (iPrint >= 99) then
  write(Label,'(A)') ' In CCrtCmp: Axyz '
  call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'R')
  call CRecPrt(Label,' ',Axyz,nZeta*3,nHer*(na+1),'I')
end if

return

end subroutine CCrtCmp
