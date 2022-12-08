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

subroutine CAssmbl(Rnxyz,Axyz,la,Bxyz,lb,nZeta,HerW,nHer)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the multipole moment *
!         matrix within the framework of the Gauss-Hermite quadrature. *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: cZero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: la, lb, nZeta, nHer
complex(kind=wp), intent(out) :: Rnxyz(nZeta*3,0:la,0:lb)
complex(kind=wp), intent(in) :: Axyz(nZeta*3,nHer,0:la), Bxyz(nZeta*3,nHer,0:lb)
real(kind=wp), intent(in) :: HerW(nHer)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iHer, iPrint, iRout
character(len=80) :: Label

iRout = 123
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In CAssmbl:HerW',' ',HerW,1,nHer)
  call CRecPrt(' In CAssmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1),'R')
  call CRecPrt(' In CAssmbl:Axyz',' ',Axyz,nZeta*3,nHer*(la+1),'I')
  call CRecPrt(' In CAssmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1),'R')
  call CRecPrt(' In CAssmbl:Bxyz',' ',Bxyz,nZeta*3,nHer*(lb+1),'I')
end if

! Initialize to zero

Rnxyz(:,:,:) = cZero

do ia=0,la
  do ib=0,lb

    ! Generate the cartesian components of the multipole moment
    ! matrix as a sum of the value of the integrand, evaluated
    ! at a root, times a weight.

    do iHer=1,nHer
      Rnxyz(:,ia,ib) = Rnxyz(:,ia,ib)+Axyz(:,iHer,ia)*Bxyz(:,iHer,ib)*HerW(iHer)
    end do

    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In CAssmbl: Rnxyz(',ia,',',ib,')'
      call CRecPrt(Label,' ',Rnxyz(:,ia,ib),nZeta,3,'R')
      call CRecPrt(Label,' ',Rnxyz(:,ia,ib),nZeta,3,'I')
    end if
  end do
end do

return

end subroutine CAssmbl
