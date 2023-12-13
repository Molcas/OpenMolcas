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

subroutine Get_maxDG(SDG,nnSkal,MxBasSh)
!***********************************************************************
!     Compute Sqrt(Abs( (mu,nu|mu,nu) ) )                              *
!     Make a list of the largest such element for each shell-pair      *
!     Store in SDG.                                                    *
!***********************************************************************

use Index_Functions, only: iTri
use Cholesky, only: iiBstR, iRS2F, iSOShl, MxOrSh, nnBstR, nnBstRT, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nnSkal
real(kind=wp), intent(inout) :: SDG(nnSkal)
integer(kind=iwp), intent(out) :: MxBasSh
integer(kind=iwp) :: iabSh, iag, iaSh, ibg, ibSh, iLoc, jRab, jSym, kRab
real(kind=wp), allocatable :: Diag(:)

!                                                                      *
!***********************************************************************
!                                                                      *
SDG(:) = Zero

iLoc = 1 ! point to 1st reduced set in index arrays
call mma_allocate(Diag,NNBSTRT(iLoc),Label='Diag')

! Read the diagonal of the integrals, (mu,nu|mu,nu)

call CHO_IODIAG(DIAG,2)

do jSym=1,nSym

  do jRab=1,nnBstR(jSym,iLoc)

    kRab = iiBstr(jSym,iLoc)+jRab ! already in 1st red set

    iag = iRS2F(1,kRab)  !global address
    ibg = iRS2F(2,kRab)

    iaSh = iSOShl(iag) ! shell to which it belongs
    ibSh = iSOShl(ibg)

    iabSh = iTri(iaSh,ibSh)

    SDG(iabSh) = max(SDG(iabSh),sqrt(abs(Diag(kRab))))

  end do  ! jRab loop
end do

call mma_deallocate(Diag)

MxBasSh = MxOrSh

return

end subroutine Get_maxDG
