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

subroutine SetUp_iSD()

use setup, only: mSkal, MxPrm
use iSD_data, only: iSD, nSD, nSkal_iSD
use k2_arrays, only: MxDij, MxFT
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iBas, iCmp, iPrim, iS, nSkal

!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(iSD)) call mma_deallocate(iSD)
call Nr_Shells(nSkal)
mSkal = nSkal
nSkal_iSD = nSkal+4  ! Add four slots for future use.
call mma_allocate(iSD,[0,nSD],[1,nSkal_iSD],label='iSD')
call Def_Shells(iSD,nSD,nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the size of and allocate auxiliary memory

MxPrm = 0
MxFT = 0
MxDij = 0
do iS=1,nSkal
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  MxPrm = max(MxPrm,iPrim)
  if (nIrrep == 1) then
    MxFT = 1 ! Dummy assignment
    MxDij = max(MxDij,iCmp**2+iPrim**2+1)
  else
    MxFT = max(MxFT,6*(iBas*iCmp)**2)
    MxDij = max(MxDij,(iBas**2+1)*iCmp**2+iPrim**2+1)
  end if
end do
MxDij = 6*nIrrep*MxDij
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine SetUp_iSD
