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

subroutine Mk_iSO2Ind(iSO2Sh,iSO2Ind,nSO,nShell)

use Basis_Info, only: nBas_Aux
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSO, iSO2Sh(nSO), nShell
integer(kind=iwp), intent(out) :: iSO2Ind(nSO)
integer(kind=iwp) :: iB, iIrrep, Ind, iSh, iSO
integer(kind=iwp), allocatable :: nTemp(:)

call mma_allocate(nTemp,nShell,Label='nTemp')

iSO = 0
do iIrrep=0,nIrrep-1

  nTemp(:) = 0
  do iB=1,nBas_Aux(iIrrep)
    iSO = iSO+1
    iSh = iSO2Sh(iSO)
    nTemp(iSh) = nTemp(iSh)+1
    Ind = nTemp(iSh)
    !write(u6,*) 'iSO,iSh,Ind=',iSO,iSh,Ind
    iSO2Ind(iSO) = Ind
  end do

end do
!call iVcPrt('iSO2Ind','(10I5)',iSO2Ind,nSO)

call mma_deallocate(nTemp)

return

end subroutine Mk_iSO2Ind
