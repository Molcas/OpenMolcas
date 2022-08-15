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

subroutine VelInt(Vxyz,Sxyz,na,nb,Beta,nZeta)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the velocity inte-   *
!         grals from the cartesian components of the overlap integals. *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: na, nb, nZeta
real(kind=wp), intent(out) :: Vxyz(nZeta,3,0:na,0:nb)
real(kind=wp), intent(in) :: Sxyz(nZeta,3,0:na,0:nb+1), Beta(nZeta)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iPrint, iRout
character(len=80) :: Label

iRout = 160
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In VelInt: Beta ',' ',Beta,nZeta,1)
end if
do ia=0,na
  do ib=0,nb
    if (ib == 0) then
      do iCar=1,3
        Vxyz(:,iCar,ia,ib) = -Beta*Two*Sxyz(:,iCar,ia,ib+1)
      end do
    else
      do iCar=1,3
        Vxyz(:,iCar,ia,ib) = real(ib,kind=wp)*Sxyz(:,iCar,ia,ib-1)-Beta*Two*Sxyz(:,iCar,ia,ib+1)
      end do
    end if

    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In VelInt: Vxyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Vxyz(:,:,ia,ib),nZeta,3)
    end if
  end do
end do

return

end subroutine VelInt
