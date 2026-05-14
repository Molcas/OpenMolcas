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
! Copyright (C) 1995, Martin Schuetz                                   *
!               1998, Roland Lindh                                     *
!               2000-2015, Steven Vancoillie                           *
!***********************************************************************

! finalize GA
subroutine GATerminate()
! SVC: This terminates the parallel runtime after which the processes
! are no longer allowed to make ga/mpi calls. When called more than
! once, the routine does nothing. This is to support early
! termination of the parallel runtime without actually exiting.
! In such a situation, when the process eventually finishes, it
! will call this routine again, but then doing nothing. Such a use
! case is e.g. when we want to terminate slave processes and only
! continue to run the master process in serial mode.

#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs
#ifdef _GA_
use Definitions, only: MPIInt
#endif
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: FirstCall = .true.
#ifdef _GA_
integer(kind=MPIInt) :: iErr
#endif

if (FirstCall) then
  FirstCall = .false.
  if (mpp_nprocs > 1) then
    call ga_terminate()
#   ifdef _GA_
    call mpi_finalize(iErr)
#   endif
  end if
end if
#endif

end subroutine GATerminate
