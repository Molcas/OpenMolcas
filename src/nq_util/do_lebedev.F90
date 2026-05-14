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

subroutine Do_Lebedev(L_Eff,mPt,R,Sym)
!***********************************************************************
!                                                                      *
!     Computes data useful for the angular quadrature.                 *
!                                                                      *
!***********************************************************************

use nq_Grid, only: Pax
use Lebedev_quadrature, only: available_table, ld_by_rule, order_table
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Four, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: L_Eff, Sym
integer(kind=iwp), intent(out) :: mPt
real(kind=wp), allocatable, intent(out) :: R(:,:)
integer(kind=iwp) :: iSet
real(kind=wp), allocatable :: TempR(:,:), TempW(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate angular grid a la Lebedev

iSet = (L_Eff-1)/2
if (available_table(iSet) /= 1) then
  write(u6,'(A,I3)') 'Failed to find a Lebedev grid of order',L_EFF
  call Abend()
end if
if ((Sym < 0) .or. (Sym > 4)) then
  write(u6,'(A,I3)') 'Symmetry of a Lebedev grid must be 0, 1, 2, 3 or 4'
  call Abend()
end if

mPt = order_table(Sym,iSet)
call mma_allocate(R,4,mPt,Label='R')
call mma_allocate(TempR,mPt,3,Label='TempR')
call mma_allocate(TempW,mPt,Label='TempW')

call ld_by_rule(Sym,iSet,TempR(:,1),TempR(:,2),TempR(:,3),TempW)
call DGEMM_('N','T',3,mPt,3,One,Pax,3,TempR,mPt,Zero,R,4)
R(4,:) = Four*Pi*TempW

call mma_deallocate(TempW)
call mma_deallocate(TempR)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Do_Lebedev
