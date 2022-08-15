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

subroutine Kntc_GIAO(Txyz,Rxyz,Wxyz,na,nb,Alpha,Beta,nZeta)
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
real(kind=wp), intent(in) :: Rxyz(nZeta,3,0:na+1,0:nb+1,0:1), Alpha(nZeta), Beta(nZeta)
real(kind=wp), intent(out) :: Txyz(nZeta,3,0:na,0:nb,0:1), Wxyz(nZeta,3,0:na,0:nb,2)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iPrint, iRout
character(len=80) :: Label

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 115
iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 99) then
  call RecPrt(' In Kntc: Alpha',' ',Alpha,nZeta,1)
  call RecPrt(' In Kntc: Beta ',' ',Beta,nZeta,1)
  do ia=0,na+1
    do ib=0,nb+1
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Rxyz(',ia,',',ib,',0)'
      call RecPrt(Label,' ',Rxyz(:,:,ia,ib,0),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Rxyz(',ia,',',ib,',1)'
      call RecPrt(Label,' ',Rxyz(:,:,ia,ib,1),nZeta,3)
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
do ia=0,na
  do ib=0,nb
    if ((ia /= 0) .and. (ib /= 0)) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib,0) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,0)+Half*ia*ib*Rxyz(:,iCar,ia-1,ib-1,0)- &
                               Beta*ia*Rxyz(:,iCar,ia-1,ib+1,0)-Alpha*ib*Rxyz(:,iCar,ia+1,ib-1,0)
        Txyz(:,iCar,ia,ib,1) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,1)+Half*ia*ib*Rxyz(:,iCar,ia-1,ib-1,1)- &
                               Beta*ia*Rxyz(:,iCar,ia-1,ib+1,1)-Alpha*ib*Rxyz(:,iCar,ia+1,ib-1,1)
        Wxyz(:,iCar,ia,ib,1) = -Two*Alpha*Rxyz(:,iCar,ia+1,ib,0)+ia*Rxyz(:,iCar,ia-1,ib,0)
        Wxyz(:,iCar,ia,ib,2) = -Two*Beta*Rxyz(:,iCar,ia,ib+1,0)+ib*Rxyz(:,iCar,ia,ib-1,0)
      end do
    else if (ia /= 0) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib,0) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,0)-Beta*ia*Rxyz(:,iCar,ia-1,ib+1,0)
        Txyz(:,iCar,ia,ib,1) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,1)-Beta*ia*Rxyz(:,iCar,ia-1,ib+1,1)
        Wxyz(:,iCar,ia,ib,1) = -Two*Alpha*Rxyz(:,iCar,ia+1,ib,0)+ia*Rxyz(:,iCar,ia-1,ib,0)
        Wxyz(:,iCar,ia,ib,2) = -Two*Beta*Rxyz(:,iCar,ia,ib+1,0)
      end do
    else if (ib /= 0) then
      do iCar=1,3
        Txyz(:,iCar,ia,ib,0) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,0)-Alpha*ib*Rxyz(:,iCar,ia+1,ib-1,0)
        Txyz(:,iCar,ia,ib,1) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,1)-Alpha*ib*Rxyz(:,iCar,ia+1,ib-1,1)
        Wxyz(:,iCar,ia,ib,1) = -Two*Alpha*Rxyz(:,iCar,ia+1,ib,0)
        Wxyz(:,iCar,ia,ib,2) = -Two*Beta*Rxyz(:,iCar,ia,ib+1,0)+ib*Rxyz(:,iCar,ia,ib-1,0)
      end do
    else
      do iCar=1,3
        Txyz(:,iCar,ia,ib,0) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,0)
        Txyz(:,iCar,ia,ib,1) = Two*Alpha*Beta*Rxyz(:,iCar,ia+1,ib+1,1)
        Wxyz(:,iCar,ia,ib,1) = -Two*Alpha*Rxyz(:,iCar,ia+1,ib,0)
        Wxyz(:,iCar,ia,ib,2) = -Two*Beta*Rxyz(:,iCar,ia,ib+1,0)
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',ib,',0)'
      call RecPrt(Label,' ',Txyz(:,:,ia,ib,0),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',ib,',1)'
      call RecPrt(Label,' ',Txyz(:,:,ia,ib,1),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Wxyz(',ia,',',ib,',1)'
      call RecPrt(Label,' ',Wxyz(:,:,ia,ib,1),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Wxyz(',ia,',',ib,',2)'
      call RecPrt(Label,' ',Wxyz(:,:,ia,ib,2),nZeta,3)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Kntc_GIAO
