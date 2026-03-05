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

! initialize GA, returns rank of process and # of processes
subroutine GAInit()
!***********************************************************************
!     purpose: initialize DGA and set the global rank and number of    *
!              processes in mpp_procid and mpp_nprocs. Then also set   *
!              the (initial) local myRank and nProcs variables.        *
!***********************************************************************

use Para_Info, only: MyRank, nProcs
#ifdef _MOLCAS_MPP_
use Para_Info, only: mpp_nprocs, mpp_procid, mpp_workshare
use Definitions, only: iwp
#endif

implicit none
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: molcas_nprocs, iRC
character(len=8) :: molcas_nprocs_env
#include "global.fh"

! SVC: bypass MPI initialization if only 1 process, this is needed for a
! specific version of GEO (so that the serial tasks which are run by MPI
! do not try to re-initialize MPI). This has the consequence that for any
! calculation where the number of processes is 1, calls to MPI/GA will fail
! at runtime (even though it will compile when inside _MOLCAS_MPP_!)
call getenvf('MOLCAS_NPROCS',molcas_nprocs_env)
if (molcas_nprocs_env(1:1) == ' ') then
  molcas_nprocs = -1
else
  read(molcas_nprocs_env,*) molcas_nprocs
end if
if (molcas_nprocs /= 1) then
# ifdef _GA_
  call mpi_init(iRC)
  call ga_initialize()
  call ga_replace_ma()
# else
  call ga_initialize()
# endif
  mpp_procid = ga_nodeid()
  mpp_nprocs = ga_nnodes()
  mpp_workshare = .true.
  ! make each slave process go to its proper work directory
  call slaveschdir(mpp_procid,iRC)
  if (iRC /= 0) call Abend()
else
  mpp_procid = 0
  mpp_nprocs = 1
  mpp_workshare = .false.
end if
MyRank = mpp_procid
nProcs = mpp_nprocs
#else
MyRank = 0
nProcs = 1
#endif

end subroutine GAInit
