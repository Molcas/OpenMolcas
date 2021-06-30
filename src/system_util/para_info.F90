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
! Copyright (C) 2014, Steven Vancoillie                                *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
!
! Common information regarding the parallel runtime
!
! GLOBAL rank and number of processes. This reflects the _entire_ MPI runtime,
! regardless of whether we are running in fake parallel or not!
!
!     mpp_procid
!     mpp_nprocs
!
! rank number of the root process, i.e., the mpp_procid of the "king"
!
!     mpp_rootid
!
! variable that describes if the MPI processes are sharing work, so if this is
! false, it means we are using a fake parallel environment and work is replicated.
!
!     mpp_workshare
!
! global variables describing the parallel environment
!
! S. Vancoillie, Aug 2014, in an effort to clean up use of parallel
! header files and global variables with overlapping functionality.
!
!     MyRank   - rank of process, starting from 0
!     nProcs   - number of processes
!
! NOTE: MyRank and nProcs are the LOCAL view of the rank and number
! of processes. That means that for a "fake" parallel run (where each
! process actually runs in serial), MyRank=0 and nProcs=1.
! This also means that in a FAKE parallel run, MyRank=0 does NOT give
! you a single "master" process or "king" of all processes!!! For that,
! you need to use the logical function KING().

module Para_Info

use Definitions, only: iwp

implicit none
private

#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: mpp_procid, mpp_nprocs
logical(kind=iwp) :: mpp_workshare
public :: mpp_procid, mpp_nprocs, mpp_workshare
#endif
integer(kind=iwp) :: MyRank, nProcs
integer(kind=iwp), parameter :: mpp_rootid = 0
public :: MyRank, nProcs, mpp_rootid, mpp_id, Is_Real_Par, King, Set_Do_Parallel

contains

function Is_Real_Par()
! SVC: determine if multiple processes are working together
  logical(kind=iwp) :: Is_Real_Par
# ifdef _MOLCAS_MPP_
  Is_Real_Par = mpp_workshare
# else
  Is_Real_Par = .false.
# endif
end function Is_Real_Par

function King()
! SVC: determine if this is the absolute master process,
!       regardless of the parallel environment
  logical(kind=iwp) :: King
# ifdef _MOLCAS_MPP_
  King = mpp_procid == mpp_rootid
# else
  King = .true.
# endif
end function King

function mpp_id()
! returns the absolute id of the process,
! regardless of the parallel environment
  integer(kind=iwp) :: mpp_id
# ifdef _MOLCAS_MPP_
  mpp_id = mpp_procid
# else
  mpp_id = mpp_rootid
# endif
end function mpp_id

subroutine Set_Do_Parallel(Par_Status)
  logical(kind=iwp), intent(in) :: Par_Status
# ifdef _MOLCAS_MPP_
  mpp_workshare = (Par_Status .and. (mpp_nprocs > 1))
  if (mpp_workshare) then
    MyRank = mpp_procid
    nProcs = mpp_nprocs
  else
    MyRank = 0
    nProcs = 1
  end if
# else
# include "macros.fh"
  unused_var(Par_Status)
  return
# endif
end subroutine Set_Do_Parallel

end module Para_Info
