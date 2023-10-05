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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_TrcIdl_Update(IAmIdle)
!
! Thomas Bondo Pedersen, May 2010.
!
! Update array for tracing idle processors

use Para_Info, only: MyRank
use Cholesky, only: Cho_Real_Par, Idle
#ifdef _DEBUGPRINT_
use Cholesky, only: LuPri, Trace_Idle
#endif
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: IAmIdle

#ifdef _DEBUGPRINT_
if ((.not. allocated(Idle)) .or. (.not. Trace_Idle)) then
  write(LuPri,'(A)') 'Cho_TrcIdl_Update should not be called in this run!'
  write(LuPri,*) 'Trace_Idle=',Trace_Idle
  call Cho_Quit('Illegal call to Cho_TrcIdl_Update',103)
end if
#endif

if (IAmIdle) then
  if (Cho_Real_Par) then
    Idle(1+myRank) = Idle(1+myRank)+1
  else
    Idle(1) = Idle(1)+1
  end if
end if

end subroutine Cho_TrcIdl_Update
