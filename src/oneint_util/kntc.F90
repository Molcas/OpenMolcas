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

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 Txyz(nZeta,3,0:na,0:nb), Sxyz(nZeta,3,0:na+1,0:nb+1), Alpha(nZeta), Beta(nZeta)
character*80 Label

iRout = 115
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In Kntc: Alpha',' ',Alpha,nZeta,1)
  call RecPrt(' In Kntc: Beta ',' ',Beta,nZeta,1)
  do ia=0,na+1
    do ib=0,nb+1
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Sxyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Sxyz(1,1,ia,ib),nZeta,3)
    end do
  end do
end if
do ia=0,na
  do ib=0,nb
    if ((ia == 0) .and. (ib == 0)) then
      do iCar=1,3
        do iZeta=1,nZeta
          Txyz(iZeta,iCar,ia,ib) = Two*Alpha(iZeta)*Beta(iZeta)*Sxyz(iZeta,iCar,ia+1,ib+1)
        end do
      end do
    else if (ia == 0) then
      do iCar=1,3
        do iZeta=1,nZeta
          Txyz(iZeta,iCar,ia,ib) = Two*Alpha(iZeta)*Beta(iZeta)*Sxyz(iZeta,iCar,ia+1,ib+1)- &
                                   Alpha(iZeta)*dble(ib)*Sxyz(iZeta,iCar,ia+1,ib-1)
        end do
      end do
    else if (ib == 0) then
      do iCar=1,3
        do iZeta=1,nZeta
          Txyz(iZeta,iCar,ia,ib) = Two*Alpha(iZeta)*Beta(iZeta)*Sxyz(iZeta,iCar,ia+1,ib+1)- &
                                   Beta(iZeta)*dble(ia)*Sxyz(iZeta,iCar,ia-1,ib+1)
        end do
      end do
    else
      do iCar=1,3
        do iZeta=1,nZeta
          Txyz(iZeta,iCar,ia,ib) = Half*dble(ia*ib)*Sxyz(iZeta,iCar,ia-1,ib-1)-Beta(iZeta)*dble(ia)*Sxyz(iZeta,iCar,ia-1,ib+1)- &
                                   Alpha(iZeta)*dble(ib)*Sxyz(iZeta,iCar,ia+1,ib-1)+ &
                                   Two*Alpha(iZeta)*Beta(iZeta)*Sxyz(iZeta,iCar,ia+1,ib+1)
        end do
      end do
    end if
    if (iPrint >= 99) then
      write(Label,'(A,I2,A,I2,A)') ' In Kntc: Txyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Txyz(1,1,ia,ib),nZeta,3)
    end if
  end do
end do

return

end subroutine Kntc
