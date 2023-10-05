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

subroutine Reord_Vk(ip_V_k,nProcs,myProc,nV_k,nV_t,nA,jSym,Array)

use Cholesky, only: InfVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nProcs, ip_V_k(nProcs), myProc, jSym, nV_k(jSym), nV_t(jSym), nA(jSym)
real(kind=wp), intent(inout) :: Array(*)
integer(kind=iwp) :: ifr, ik, iSym, ito, jOff, kOff, nAV_t
real(kind=wp), allocatable :: Scr(:)

nAV_t = 0
do iSym=1,jSym
  nAV_t = nAV_t+nA(iSym)*nV_t(iSym)
end do
call mma_allocate(Scr,nAV_t,Label='Scr')
Scr(:) = Zero

! On input Array first blocked over the processes
!    pointer to the block is ip_V_K(i)
!    Each block is symmetry blocked
!         Each symmetry is nA(iSym)*nV_k(iSym)
!
! Scr is also symmetry blocked
!    Each symmetry is nA(iSym)*nV_t(iSym)
!
! InfVec(ik,5,iSym) translates the local ik'th index into
! the global index of V

jOff = 0
kOff = 0
do iSym=1,jSym
  do ik=1,nV_k(iSym)   ! loop over the local vector

    ifr = ip_V_k(myProc)+jOff+nA(iSym)*(ik-1)
    ito = 1+kOff+nA(iSym)*(InfVec(ik,5,iSym)-1)
    Scr(ito:ito+nA(iSym)-1) = Array(ifr:ifr+nA(iSym)-1)

  end do
  jOff = jOff+nA(iSym)*nV_k(iSym)
  kOff = kOff+nA(iSym)*nV_t(iSym)
end do

! Copy Scr => Array
Array(ip_V_k(1):ip_V_k(1)+nAV_t-1) = Scr(1:nAV_t)

! Make a global add to get the contributions from all nodes.

call GADGOP(Array(ip_V_k(1)),nAV_t,'+')

call mma_deallocate(Scr)

return

end subroutine Reord_Vk
