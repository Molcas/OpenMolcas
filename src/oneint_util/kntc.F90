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

subroutine Kntc(Txyz,Sxyz,na,nb,Alpha,Beta,nZeta)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the kinetic energy   *
!         integral from the cartesian components of the overlap        *
!         integral.                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: na, nb, nZeta
real(kind=wp), intent(out) :: Txyz(nZeta,3,0:na,0:nb)
real(kind=wp), intent(in) :: Sxyz(nZeta,3,0:na+1,0:nb+1), Alpha(nZeta), Beta(nZeta)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iPrint, iRout
character(len=80) :: Label

iRout = 115
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In Kntc: Alpha',' ',Alpha,nZeta,1)
  call RecPrt(' In Kntc: Beta ',' ',Beta,nZeta,1)
  do ia=0,na+1
    do ib=0,nb+1
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Sxyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Sxyz(:,:,ia,ib),nZeta,3)
    end do
  end do
end if
do ia=0,na
  do ib=0,nb
    if ((ia /= 0) .and. (ib /= 0)) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib) = Two*Alpha*Beta*Sxyz(:,iCar,ia+1,ib+1)+Half*real(ia*ib,kind=wp)*Sxyz(:,iCar,ia-1,ib-1)- &
                             Alpha*real(ib,kind=wp)*Sxyz(:,iCar,ia+1,ib-1)-Beta*real(ia,kind=wp)*Sxyz(:,iCar,ia-1,ib+1)
      end do
    else if (ia /= 0) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib) = Two*Alpha*Beta*Sxyz(:,iCar,ia+1,ib+1)-Beta*real(ia,kind=wp)*Sxyz(:,iCar,ia-1,ib+1)
      end do
    else if (ib /= 0) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib) = Two*Alpha*Beta*Sxyz(:,iCar,ia+1,ib+1)-Alpha*real(ib,kind=wp)*Sxyz(:,iCar,ia+1,ib-1)
      end do
    else
      do iCar=1,3
        Txyz(:,iCar,ia,ib) = Two*Alpha*Beta*Sxyz(:,iCar,ia+1,ib+1)
      end do
    end if
    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Txyz(:,:,ia,ib),nZeta,3)
    end if
  end do
end do

return

end subroutine Kntc
