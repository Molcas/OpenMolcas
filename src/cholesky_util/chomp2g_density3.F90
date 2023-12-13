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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_Density3(irc,CMO)
!***********************************************************************
!     Jonas Bostrom, March 2010.                                       *
!                                                                      *
!     Purpose: Finalize MP2 Density.                                   *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Cholesky, only: nSym
use Cholesky_procedures, only: Build_Mp2Dens
use ChoMP2, only: MP2D, MP2D_e, MP2W, MP2W_e, nDel, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp) :: i, iSym, j, lTriDens, nOrbAll(8)
real(kind=wp), allocatable :: AOTriDens(:), WAOTriDens(:)

!                                                                      *
!***********************************************************************
!                                                                      *
irc = 0

nOrbAll(:) = nOrb(:)+nDel(:)

lTriDens = 0
do iSym=1,nSym
  lTriDens = lTriDens+nTri_Elem(nOrbAll(iSym))
end do

do iSym=1,nSym
  do i=1,nOrbAll(iSym)
    do j=1,nOrbAll(iSym)
      if ((i <= nOrb(iSym)) .and. (j <= nOrb(iSym))) then
        MP2D_e(iSym)%A(i,j) = MP2D(iSym)%A(i,j)
        MP2W_e(iSym)%A(i,j) = MP2W(iSym)%A(i,j)
      else
        MP2D_e(iSym)%A(i,j) = Zero
        MP2W_e(iSym)%A(i,j) = Zero
      end if
    end do
  end do
end do

call mma_allocate(AOTriDens,lTriDens,Label=' AOTriDens')
call mma_allocate(WAOTriDens,lTriDens,Label='WAOTriDens')
AOTriDens(:) = Zero
WAOTriDens(:) = Zero

call Build_Mp2Dens(AOTriDens,lTriDens,MP2D_e,CMO,nSym,nOrbAll,.true.)
call Build_Mp2Dens(WAOTriDens,lTriDens,MP2W_e,CMO,nSym,nOrbAll,.false.)

call Put_dArray('D1aoVar',AOTriDens,lTriDens)
call Put_dArray('FockOcc',WAOTriDens,lTriDens)

call mma_deallocate(AOTriDens)
call mma_deallocate(WAOTriDens)

end subroutine ChoMP2g_Density3
