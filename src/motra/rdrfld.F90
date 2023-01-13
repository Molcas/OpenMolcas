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

subroutine RdRfld()

use motra_global, only: HOne, nBas, nSym, PotNuc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iSym, nTemp
real(kind=wp) :: ERFself
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Temp(:)

!----------------------------------------------------------------------*
! If this is a perturbative reaction field calculation then            *
! modifiy the one-electron Hamiltonian by the reaction field and       *
! the nuclear attraction by the cavity self-energy                     *
!----------------------------------------------------------------------*
nTemp = 0
do iSym=1,nSym
  nTemp = nTemp+nBas(iSym)*(nBas(iSym)+1)/2
end do
call mma_allocate(Temp,nTemp,label='RFFLD')
call f_Inquire('RUNOLD',Found)
if (Found) call NameRun('RUNOLD')
call get_dscalar('RF Self Energy',ERFself)
PotNuc = PotNuc+ERFself
call get_darray('Reaction field',Temp,nTemp)
if (Found) call NameRun('#Pop')
call Daxpy_(nTemp,One,Temp,1,HOne,1)
call mma_deallocate(Temp)
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
return

end subroutine RdRfld
