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

subroutine INTCTL2(CMO,nCMO,DREF,nDREF,FIFA,NFIFA,FIMO,nFIMO)

use PrintLevel, only: DEBUG
use caspt2_global, only: do_grad, FIFA_all, FIMO_all, iPrGlb, nStpGrd
use caspt2_module, only: nBTri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCMO, nDREF, NFIFA, nFIMO
real(kind=wp), intent(in) :: CMO(nCMO), DREF(nDREF)
real(kind=wp), intent(out) :: FIFA(NFIFA), FIMO(nFIMO)
real(kind=wp), allocatable :: FFIAO(:), FAAO(:)
logical(kind=iwp), parameter :: IF_TRNSF = .false.

! Compute using Cholesky vectors.
! Frozen + inactive and active Fock matrix in AO basis:
! Changed to calculate frozen and inactive Fock separately to frozen + inactive, because there are no reasons that we do so
call mma_allocate(FFIAO,NBTRI,LABEL='FFIAO')
call mma_allocate(FAAO,NBTRI,LABEL='FAAO')

! tracho2 makes many allocations but should deallocate everything
! before its return.
if (IPRGLB >= DEBUG) write(u6,*) ' INTCTL2 calling TRACHO2...'

call TraCho2(CMO,nCMO,DREF,nDREF,FFIAO,FAAO,IF_TRNSF)

if (IPRGLB >= DEBUG) write(u6,*) ' INTCTL2 back from TRACHO2.'
! All extra allocations inside tracho2 should now be gone.

! For gradient calculation, it is good to have FFIAO and FAAO
if (do_grad .or. (nStpGrd == 2)) then

  !! FFIAO has one-electron Hamiltonian
  FIMO_all(1:NBTri) = FFIAO(1:NBTri)
  FIFA_all(1:NBTri) = FIMO_all(1:NBTri)+FAAO(1:NBTri)

end if

! Transform to MO basis: generating FIMO and FIFA.
call FMat_Cho(CMO,nCMO,FFIAO,FAAO,FIMO,nFIMO,FIFA,nFIFA)

call mma_deallocate(FFIAO)
call mma_deallocate(FAAO)

end subroutine INTCTL2
