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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module EDA_Parallel

use Definitions, only: iwp, u6
use Para_Info, only: nProcs, MyRank, Is_Real_Par

implicit none
private

public :: EDA_Sync_Error

contains

subroutine EDA_Sync_Error(Code,Message)

  integer(kind=iwp), intent(inout) :: Code
  character(len=*), intent(inout) :: Message
  logical :: FailedHere
  integer(kind=iwp) :: AnyError

  FailedHere = (Code /= 0)
  AnyError = 0
  if (FailedHere) AnyError = 1
  if ((nProcs > 1) .and. Is_Real_Par()) call gaIgOP_SCAL(AnyError,'max')
  if (AnyError /= 0) then
    if (FailedHere) write(u6,'(A,I6,2A)') ' EDA rank ',MyRank,': ',trim(Message)
    Code = 1
    if (.not. FailedHere) Message = 'EDA failed on another MPI rank: inspect rank diagnostics.'
  end if

end subroutine EDA_Sync_Error

end module EDA_Parallel
