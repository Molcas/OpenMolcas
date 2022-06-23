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

subroutine CVelInt(Vxyz,Sxyz,na,nb,Alpha,Beta,nZeta)
!***********************************************************************
!                                                                      *
! Object: to assemble the cartesian components of the velocity inte-   *
!         grals from the cartesian components of the overlap integals. *
!                                                                      *
! Called from: PrpInt                                                  *
!                                                                      *
! Calling    : CRecPrt                                                 *
!               RecPrt                                                 *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: na, nb, nZeta
complex(kind=wp) :: Vxyz(nZeta,3,0:na,0:nb,2), Sxyz(nZeta,3,0:na+1,0:nb+1)
real(kind=wp) :: Alpha(nZeta), Beta(nZeta)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iPrint, iRout, iZeta
character(len=80) :: Label

iRout = 160
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In CVelInt: Beta ',' ',Beta,nZeta,1)
end if
do ia=0,na
  do ib=0,nb
    if ((ia /= 0) .and. (ib /= 0)) then
      do iCar=1,3
        do iZeta=1,nZeta
          Vxyz(iZeta,iCar,ia,ib,1) = real(ia,kind=wp)*Sxyz(iZeta,iCar,ia-1,ib)-Alpha(iZeta)*Two*Sxyz(iZeta,iCar,ia+1,ib)
          Vxyz(iZeta,iCar,ia,ib,2) = real(ib,kind=wp)*Sxyz(iZeta,iCar,ia,ib-1)-Beta(iZeta)*Two*Sxyz(iZeta,iCar,ia,ib+1)
        end do
      end do
    else if ((ia == 0) .and. (ib /= 0)) then
      do iCar=1,3
        do iZeta=1,nZeta
          Vxyz(iZeta,iCar,ia,ib,1) = -Alpha(iZeta)*Two*Sxyz(iZeta,iCar,ia+1,ib)
          Vxyz(iZeta,iCar,ia,ib,2) = real(ib,kind=wp)*Sxyz(iZeta,iCar,ia,ib-1)-Beta(iZeta)*Two*Sxyz(iZeta,iCar,ia,ib+1)
        end do
      end do
    else if ((ia /= 0) .and. (ib == 0)) then
      do iCar=1,3
        do iZeta=1,nZeta
          Vxyz(iZeta,iCar,ia,ib,1) = real(ia,kind=wp)*Sxyz(iZeta,iCar,ia-1,ib)-Alpha(iZeta)*Two*Sxyz(iZeta,iCar,ia+1,ib)
          Vxyz(iZeta,iCar,ia,ib,2) = -Beta(iZeta)*Two*Sxyz(iZeta,iCar,ia,ib+1)
        end do
      end do
    else
      do iCar=1,3
        do iZeta=1,nZeta
          Vxyz(iZeta,iCar,ia,ib,1) = -Alpha(iZeta)*Two*Sxyz(iZeta,iCar,ia+1,ib)
          Vxyz(iZeta,iCar,ia,ib,2) = -Beta(iZeta)*Two*Sxyz(iZeta,iCar,ia,ib+1)
        end do
      end do
    end if

    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In CVelInt: Vxyz(',ia,',',ib,',1)'
      call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,1),nZeta,3,'R')
      call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,1),nZeta,3,'I')
      write(Label,'(A,I2,A,I2,A)') ' In CVelInt: Vxyz(',ia,',',ib,',2)'
      call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,2),nZeta,3,'R')
      call CRecPrt(Label,' ',Vxyz(1,1,ia,ib,2),nZeta,3,'I')
    end if
  end do
end do

return

end subroutine CVelInt
