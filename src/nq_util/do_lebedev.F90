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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Do_Lebedev(L_Eff,mPt,R)
!***********************************************************************
!                                                                      *
!     Computes data useful for the angular quadrature.                 *
!                                                                      *
!***********************************************************************

use nq_Grid, only: Pax
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: L_Eff
integer(kind=iwp), intent(out) :: mPt
real(kind=wp), allocatable, intent(out) :: R(:,:)
integer(kind=iwp) :: iSet, nPt
integer(kind=iwp), parameter :: nSet = 11, Lebedev_order(nSet) = [5,7,11,17,23,29,35,41,47,53,59], &
                                Lebedev_npoints(nSet) = [14,26,50,110,194,302,434,590,770,974,1202]
real(kind=wp), allocatable :: TempR(:,:), TempW(:)

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
      write(u6,*) 'nPt=',nPt
      write(u6,*) 'mPt=',mPt
      call Abend()
    end if

    call DGEMM_('N','N',3,nPt,3,One,Pax,3,TempR,3,Zero,R,4)
    R(4,:) = Four*Pi*TempW

    call mma_deallocate(TempW)
    call mma_deallocate(TempR)

    return

  end if
end do
write(u6,'(A,I3)') 'Failed to find a Lebedev grid of order',L_EFF
write(u6,'(A)') 'Available orders are:'
write(u6,'(11(1X,I3))') Lebedev_order(:)
call Abend()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Do_Lebedev
