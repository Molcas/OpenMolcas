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

! Common information regarding the parallel runtime

! GLOBAL rank and number of processes. This reflects the _entire_ MPI runtime,
! regardless of whether we are running in fake parallel or not!

!     mpp_procid
!     mpp_nprocs

! rank number of the root process, i.e., the mpp_procid of the "king"

!     mpp_rootid

! variable that describes if the MPI processes are sharing work, so if this is
! false, it means we are using a fake parallel environment and work is replicated.

!     mpp_workshare

! global variables describing the parallel environment
!
! S. Vancoillie, Aug 2014, in an effort to clean up use of parallel
! header files and global variables with overlapping functionality.

!     MyRank   - rank of process, starting from 0
!     nProcs   - number of processes

! NOTE: MyRank and nProcs are the LOCAL view of the rank and number
! of processes. That means that for a "fake" parallel run (where each
! process actually runs in serial), MyRank=0 and nProcs=1.
! This also means that in a FAKE parallel run, MyRank=0 does NOT give
! you a single "master" process or "king" of all processes!!! For that,
! you need to use the logical function KING().

Module Para_Info

Implicit None
Private
#ifdef _MOLCAS_MPP_
Integer :: mpp_procid, mpp_nprocs
Logical :: mpp_workshare
Public :: mpp_procid, mpp_nprocs, mpp_workshare
#endif
Integer :: MyRank, nProcs
Integer, Parameter :: mpp_rootid = 0
Public :: MyRank, nProcs, mpp_rootid, mpp_id, Is_Real_Par, King, Set_Do_Parallel

Contains

Logical Function Is_Real_Par()
! SVC: determine if multiple processes are working together
#ifdef _MOLCAS_MPP_
  Is_Real_Par = mpp_workshare
#else
  Is_Real_Par = .FALSE.
#endif
End Function Is_Real_Par

Logical Function King()
! SVC: determine if this is the absolute master process,
!       regardless of the parallel environment
#ifdef _MOLCAS_MPP_
  King = mpp_procid == mpp_rootid
#else
  King = .TRUE.
#endif
End Function King

Integer Function mpp_id()
! returns the absolute id of the process,
! regardless of the parallel environment
#ifdef _MOLCAS_MPP_
  mpp_id = mpp_procid
#else
  mpp_id = mpp_rootid
#endif
End Function mpp_id

Subroutine Set_Do_Parallel(Par_Status)
  ! Set to optional to avoid unused argument warnings
  Logical, Intent(In), Optional :: Par_Status
#ifdef _MOLCAS_MPP_
  If (Present(Par_Status)) mpp_workshare = (Par_Status .and. (mpp_nprocs > 1))
  If (mpp_workshare) Then
    MyRank = mpp_procid
    nProcs = mpp_nprocs
  Else
    MyRank = 0
    nProcs = 1
  End If
#else
  Return
#endif
End Subroutine Set_Do_Parallel

End Module Para_Info
