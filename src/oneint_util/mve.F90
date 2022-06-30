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
! Copyright (C) 1990,1998, Roland Lindh                                *
!               1998, Samuel Mikes                                     *
!***********************************************************************

subroutine MVe(rV2Int,rV4Int,Sxyz,na,nb,Alpha,Beta,nZeta)
!***********************************************************************
!                                                                      *
! Object: to compute intermediate integrals for the evaluation of the  *
!          mass velocity integrals.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!                                                                      *
!             Correct out of bound reference, February '98, Samuel     *
!             Mikes and Roland Lindh.                                  *
!***********************************************************************

use Constants, only: Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: na, nb, nZeta
real(kind=wp), intent(out) :: rV2Int(nZeta,3,0:na,0:nb,2), rV4Int(nZeta,3,0:na,0:nb)
real(kind=wp), intent(in) :: Sxyz(nZeta,3,0:na+2,0:nb+2), Alpha(nZeta), Beta(nZeta)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iCar, iPrint, iRout
character(len=80) :: Label

iRout = 192
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  call RecPrt(' In MVe: Alpha',' ',Alpha,nZeta,1)
  call RecPrt(' In MVe: Beta ',' ',Beta,nZeta,1)
  do ib=0,nb+2
    do ia=0,na+2
      write(Label,'(A,I2,A,I2,A)') ' In MVe: Sxyz(',ia,',',ib,')'
      call RecPrt(Label,' ',Sxyz(:,:,ia,ib),nZeta,3)
    end do
  end do

end if
do ib=0,nb
  do ia=0,na
    do iCar=1,3
      rV2Int(:,iCar,ia,ib,1) = Four*Alpha**2*Sxyz(:,iCar,ia+2,ib)-Two*Alpha*real(2*ia+1,kind=wp)*Sxyz(:,iCar,ia,ib)
      if (ia >= 2) rV2Int(:,iCar,ia,ib,1) = rV2Int(:,iCar,ia,ib,1)+real(ia*(ia-1),kind=wp)*Sxyz(:,iCar,ia-2,ib)

      rV2Int(:,iCar,ia,ib,2) = Four*Beta**2*Sxyz(:,iCar,ia,ib+2)-Two*Beta*real(2*ib+1,kind=wp)*Sxyz(:,iCar,ia,ib)
      if (ib >= 2) rV2Int(:,iCar,ia,ib,2) = rV2Int(:,iCar,ia,ib,2)+real(ib*(ib-1),kind=wp)*Sxyz(:,iCar,ia,ib-2)

      rV4Int(:,iCar,ia,ib) = Four*Alpha**2*Four*Beta**2*Sxyz(:,iCar,ia+2,ib+2)- &
                             Four*Alpha**2*Two*Beta*real(2*ib+1,kind=wp)*Sxyz(:,iCar,ia+2,ib)- &
                             Four*Beta**2*Two*Alpha*real(2*ia+1,kind=wp)*Sxyz(:,iCar,ia,ib+2)+ &
                             Two*Alpha*real(2*ia+1,kind=wp)*Two*Beta*real(2*ib+1,kind=wp)*Sxyz(:,iCar,ia,ib)
      if (ia >= 2) rV4Int(:,iCar,ia,ib) = rV4Int(:,iCar,ia,ib)+real(ia*(ia-1),kind=wp)* &
                                          (Four*Beta**2*Sxyz(:,iCar,ia-2,ib+2)-Two*Beta*real(2*ib+1,kind=wp)*Sxyz(:,iCar,ia-2,ib))
      if (ib >= 2) rV4Int(:,iCar,ia,ib) = rV4Int(:,iCar,ia,ib)+real(ib*(ib-1),kind=wp)* &
                                          (Four*Alpha**2*Sxyz(:,iCar,ia+2,ib-2)-Two*Alpha*real(2*ia+1,kind=wp)*Sxyz(:,iCar,ia,ib-2))
      if ((ia >= 2) .and. (ib >= 2)) rV4Int(:,iCar,ia,ib) = rV4Int(:,iCar,ia,ib)+ &
                                                            real(ia*(ia-1)*ib*(ib-1),kind=wp)*Sxyz(:,iCar,ia-2,ib-2)
    end do
  end do
end do

if (iPrint >= 99) then
  do ib=0,nb
    do ia=0,na
      write(Label,'(A,I2,A,I2,A)') 'In MVe: rV2Int(',ia,',',ib,',1)'
      call RecPrt(Label,' ',rV2Int(:,:,ia,ib,1),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') 'In MVe: rV2Int(',ia,',',ib,',2)'
      call RecPrt(Label,' ',rV2Int(:,:,ia,ib,2),nZeta,3)
      write(Label,'(A,I2,A,I2,A)') 'In MVe: rV4Int(',ia,',',ib,')'
      call RecPrt(Label,' ',rV4Int(:,:,ia,ib),nZeta,3)
    end do
  end do
end if

return

end subroutine MVe
