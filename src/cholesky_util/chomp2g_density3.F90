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

use ChoMP2, only: MP2D, MP2W, MP2W_e, MP2D_e
use Constants
use stdalloc
use ChoMP2g

implicit real*8(a-h,o-z)
integer irc
real*8 CMO(*)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
character(len=8), parameter :: ThisNm = 'Density3'
character(len=16), parameter :: SecNam = 'ChoMP2g_Density3'
integer nOccAll(8), nOrbAll(8)
real*8, allocatable :: AOTriDens(:), WAOTriDens(:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Build_Mp2Dens(TriDens,nTriDens,MP2X_e,CMO,mSym,nOrbAll,nOccAll,Diagonalize)

    use ChoMP2, only: Pointer_2D

    integer, intent(in) :: nTriDens
    real*8, intent(inout) :: TriDens(nTriDens)
    type(Pointer_2D), intent(in) :: MP2X_e(8)
    real*8, intent(in) :: CMO(*)
    integer, intent(in) :: mSym
    integer, intent(in) :: nOrbAll(8), nOccAll(8)
    logical, intent(in) :: Diagonalize
  end subroutine Build_Mp2Dens
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
irc = 0

do iSym=1,8
  nOccAll(iSym) = nOcc(iSym)+nFro(iSym)
  nOrbAll(iSym) = nOrb(iSym)+nDel(iSym)
end do

lTriDens = 0
do iSym=1,nSym
  lTriDens = lTriDens+nOrbAll(iSym)*(nOrbAll(iSym)+1)/2
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

call Build_Mp2Dens(AOTriDens,lTriDens,MP2D_e,CMO,nSym,nOrbAll,nOccAll,.true.)
call Build_Mp2Dens(WAOTriDens,lTriDens,MP2W_e,CMO,nSym,nOrbAll,nOccAll,.false.)

call Put_dArray('D1aoVar',AOTriDens,lTriDens)
call Put_dArray('FockOcc',WAOTriDens,lTriDens)

call mma_deallocate(AOTriDens)
call mma_deallocate(WAOTriDens)

end subroutine ChoMP2g_Density3
