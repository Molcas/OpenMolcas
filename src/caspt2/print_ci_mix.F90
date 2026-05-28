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
! Copyright (C) 1997, Per Ake Malmqvist                                *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Print_CI_Mix(EigVec)

use RefWfn, only: iadr15, refwfn_active, refwfn_close, refwfn_filename, refwfn_id, refwfn_is_h5
use caspt2_module, only: CIThr, nConf, nState, STSym
#ifdef _HDF5_
use mh5, only: mh5_fetch_dset, mh5_open_file_r
use caspt2_module, only: mState
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: EigVec(nState,nState)
integer(kind=iwp) :: iDisk, iiState, iState
#ifdef _HDF5_
integer(kind=iwp) :: jSNum
#endif
logical(kind=iwp) :: Close_refwfn
real(kind=wp), allocatable :: cCI(:), mCI(:)

call mma_allocate(mCI,nConf,Label='MixCICoeff')
call mma_allocate(cCI,nConf,Label='CICoeff')

Close_refwfn = .false.
if (.not. refwfn_active) then
  ! bypass refwfn_open, because we don't want to set global stuff
  if (refwfn_is_h5) then
#   ifdef _HDF5_
    refwfn_id = mh5_open_file_r(refwfn_filename)
#   else
    ! This should never happen
    call AbEnd()
#   endif
  else
    refwfn_id = 15
    call DAName(refwfn_id,refwfn_filename)
  end if
  Close_refwfn = .true.
end if

call CollapseOutput(1,'Mixed CI coefficients:')

write(u6,*)
write(u6,*) ' The original CI arrays are now mixed as linear'
write(u6,*) ' combinations, given by the eigenvectors.'
write(u6,*)

do iState=1,nState
  mCI(:) = Zero
  iDisk = iAdr15(4)
  do iiState=1,nState
    if (refwfn_is_h5) then
#     ifdef _HDF5_
      jSNum = mState(iiState)
      call mh5_fetch_dset(refwfn_id,'CI_VECTORS',cCI,[nConf,1],[0,jSNum-1])
#     else
      ! This should never happen
      call AbEnd()
#     endif
    else
      call dDAFile(refwfn_id,2,cCI,nConf,iDisk)
    end if
    mCI(:) = EigVec(iiState,iState)*cCI(:)
  end do
  write(u6,'(1X,A,I3)') ' The CI coefficients for the MIXED state nr. ',iState
  call PrWf_CP2(stSym,nConf,mCI,CITHR)
end do

call CollapseOutput(0,'Mixed CI coefficients:')
write(u6,*)

if (Close_refwfn) call refwfn_close()

call mma_deallocate(mCI)
call mma_deallocate(cCI)

end subroutine Print_CI_Mix
