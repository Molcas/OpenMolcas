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

subroutine Cho_TrcIdl_Report()
!
! Thomas Bondo Pedersen, May 2010.
!
! Report idle status for all processors

use Para_Info, only: nProcs
use Cholesky, only: Cho_Real_Par, Idle, LuPri
#ifdef _DEBUGPRINT_
use Cholesky, only: Trace_Idle
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, nIdle
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: l_Idle
#endif
integer(kind=iwp), allocatable :: TIloc(:)

#ifdef _DEBUGPRINT_
if ((.not. allocated(Idle)) .or. (.not. Trace_Idle)) then
  write(LuPri,'(A)') 'Cho_TrcIdl_Report should not be called in this run!'
  write(LuPri,*) 'Trace_Idle=',Trace_Idle
  call Cho_Quit('Illegal call to Cho_TrcIdl_Report',103)
end if
#endif

if (Cho_Real_Par) then
# ifdef _DEBUGPRINT_
  l_Idle = 0
  if (allocated(Idle)) l_Idle = size(Idle)
  if (l_Idle < nProcs) then
    write(LuPri,'(A)') 'Error detected in Cho_TrcIdl_Report: l_Idle < nProcs'
    write(LuPri,*) 'Trace_Idle=',Trace_Idle
    call Cho_Quit('Cho_TrcIdle_Report: l_Idle not properly set!',103)
  end if
# endif
  call mma_allocate(TIloc,[0,nProcs-1],Label='TILoc')
  TIloc(:) = Idle(1:nProcs)
  call Cho_GAIGOp(TIloc,nProcs,'+')
  nIdle = 0
  do i=0,nProcs-1
    nIdle = nIdle+min(TIloc(i),1)
  end do
  if (nIdle == 0) then
    write(LuPri,'(A)') 'No idle procs to report'
  else
    write(LuPri,'(I4,A,I4,A,F7.2,A)') nIdle,' of',nProcs,' procs have been idle (', &
                                      1.0e2_wp*real(nIdle,kind=wp)/real(nProcs,kind=wp),' %)'
    write(LuPri,'(A)') 'List of idle procs:'
    do i=0,nProcs-1
      if (TIloc(i) > 0) write(LuPri,'(I4,A,I8,A)') i,' (Idle counter:',TIloc(i),')'
    end do
  end if
  call mma_deallocate(TIloc)
else
  if (Idle(1) == 0) then
    write(LuPri,'(A)') 'No idle procs to report!'
  else
    write(LuPri,'(A,I8,A)') 'Proc 0 has been idle',Idle(1),' times'
  end if
end if
call XFlush(LuPri)

end subroutine Cho_TrcIdl_Report
