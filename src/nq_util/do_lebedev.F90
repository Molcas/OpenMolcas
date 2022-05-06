!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Do_Lebedev(L_Eff,mPt,R)
!***********************************************************************
!                                                                      *
!     Computes datas useful for the angular quadrature.                *
!                                                                      *
!***********************************************************************

use nq_Grid, only: Pax
use nq_Info

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
parameter(nSet=11)
integer Lebedev_order(nSet), Lebedev_npoints(nSet)
data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
data Lebedev_npoints/14,26,50,110,194,302,434,590,770,974,1202/
real*8, allocatable :: TempR(:,:), TempW(:)
real*8, allocatable :: R(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate angular grid a la Lebedev

do iSet=1,nSet
  if (Lebedev_order(iSet) == L_Eff) then
    mPt = Lebedev_npoints(iSet)
    call mma_allocate(R,4,mPt,Label='R')
    call mma_allocate(TempR,3,mPt,Label='TempR')
    call mma_allocate(TempW,mPt,Label='TempW')

    call Lebedev(TempR,TempW,nPt,mPt,L_Eff)
    if (nPt /= mPt) then
      call WarningMessage(2,'Lebedev_Grid: nPt /= mPt')
      write(6,*) 'nPt=',nPt
      write(6,*) 'mPt=',mPt
      call Abend()
    end if
    call DScal_(nPt,Four*Pi,TempW,1)

    call DGEMM_('N','N',3,nPt,3,1.0d0,Pax,3,TempR,3,0.0d0,R,4)
    call dcopy_(nPt,TempW,1,R(4,1),4)

    call mma_deallocate(TempW)
    call mma_deallocate(TempR)

    return

  end if
end do
write(6,'(A,I3)') 'Failed to find a Lebedev grid of order',L_EFF
write(6,'(A)') 'Available orders are:'
write(6,'(11(1X,I3))') (Lebedev_order(i),i=1,nSet)
call Abend()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Do_Lebedev
