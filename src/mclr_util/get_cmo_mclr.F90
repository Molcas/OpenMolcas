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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

subroutine Get_CMO_MCLR(CMO,Length)
! Read the MO coefficients of the last RASSCF iteration from the runfile
! as nBas x nOrb blocks (the format used in MCLR). RASSCF stores them
! as nBas x nBas blocks, including the deleted orbitals as the last columns
! of each symmetry block; these are dropped here.

use Orbital_Spaces, only: nSym
use input_mclr, only: nBas, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Length
real(kind=wp), intent(out) :: CMO(Length)
integer(kind=iwp) :: iSym, ipF, ipO, nData
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Full(:)

call Qpg_dArray('Last orbitals',Found,nData)
if (Found .and. (nData == Length)) then
  call Get_dArray('Last orbitals',CMO,Length)
else if (Found .and. (nData == sum(nBas(1:nSym)**2))) then
  call mma_allocate(Full,nData,Label='Full')
  call Get_dArray('Last orbitals',Full,nData)
  ipF = 0
  ipO = 0
  do iSym=1,nSym
    CMO(ipO+1:ipO+nBas(iSym)*nOrb(iSym)) = Full(ipF+1:ipF+nBas(iSym)*nOrb(iSym))
    ipF = ipF+nBas(iSym)**2
    ipO = ipO+nBas(iSym)*nOrb(iSym)
  end do
  call mma_deallocate(Full)
else
  write(u6,*) 'Get_CMO_MCLR: unexpected length of Last orbitals:',nData,Length
  call Abend()
end if

end subroutine Get_CMO_MCLR
