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

subroutine INTCTL2(CMO,nCMO,DREF,nDREF,FIFA,NFIFA,HONE,nHONE,FIMO,nFIMO)

use caspt2_global, only: iPrGlb
use caspt2_global, only: do_grad, nStpGrd, FIMO_all, FIFA_all
use PrintLevel, only: DEBUG
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: nBTri
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: nCMO, nDREF, NFIFA, nHONE, nFIMO
real(kind=wp), intent(in) :: CMO(nCMO), DREF(nDREF), HONE(nHONE)
real(kind=wp), intent(out) :: FIFA(NFIFA), FIMO(nFIMO)
logical(kind=iwp), parameter :: IF_TRNSF = .false.
real(kind=wp), allocatable :: FFAO(:), FIAO(:), FAAO(:)

! Compute using Cholesky vectors.
! Frozen, inactive and active Fock matrix in AO basis:
call mma_allocate(FFAO,NBTRI,LABEL='FFAO')
call mma_allocate(FIAO,NBTRI,LABEL='FIAO')
call mma_allocate(FAAO,NBTRI,LABEL='FAAO')

! tracho2 makes many allocations but should deallocate everything
! before its return.
if (IPRGLB >= DEBUG) then
  write(6,*) ' INTCTL2 calling TRACHO2...'
  call XFLUSH(6)
end if

call TraCho2(CMO,nCMO,DREF,nDREF,FFAO,FIAO,FAAO,IF_TRNSF)

if (IPRGLB >= DEBUG) then
  write(6,*) ' INTCTL2 back from TRACHO2.'
  call XFLUSH(6)
end if
! All extra allocations inside tracho2 should now be gone.

! For gradient calculation, it is good to have FIAO and FAAO
if (do_grad .or. (nStpGrd == 2)) then

  !! FFAO has one-electron Hamiltonian
  FIMO_all(1:NBTri) = FFAO(1:NBTri)+FIAO(1:NBTri)
  FIFA_all(1:NBTri) = FIMO_all(1:NBTri)+FAAO(1:NBTri)

end if

! Transform to MO basis: generating HONE, FIMO and FIFA
call FMat_Cho(CMO,nCMO,FIAO,FAAO,HONE,nHONE,FIMO,nFIMO,FIFA,nFIFA)

call mma_deallocate(FFAO)
call mma_deallocate(FIAO)
call mma_deallocate(FAAO)

end subroutine INTCTL2
