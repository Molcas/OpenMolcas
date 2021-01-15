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
! stefan: interface wrapper to the QCMaquis interface controller (qcmaquis_interface_ctl)

module qcmaquis_interface_wrapper

#ifdef _DMRG_
 use qcmaquis_interface_main
 use qcmaquis_interface_cfg

implicit none

  !> public
  public dmrg_interface_ctl

contains

  subroutine dmrg_interface_ctl(                                  &
                                task,                             &
                                x1,                               &
                                x2,                               &
                                x3,                               &
                                x4,                               &
                                energy,                           &
                                ndim,                             &
                                mdim,                             &
                                odim,                             &
                                pdim,                             &
                                state,                            &
                                stateL,                           &
                                msproj,                           &
                                msprojL,                          &
                                multiplet,                        &
                                multipletL,                       &
                                rdm1,                             &
                                rdm2,                             &
                                rdm3,                             &
                                rdm4,                             &
                                Key_CION,                         &
                                IterSCF,                          &
                                checkpoint1,                      &
                                checkpoint2                       &
                               )

      character(len=8),                intent(in)    :: task
      real*8 , optional, dimension(*), intent(inout) :: x1
      real*8 , optional, dimension(*), intent(inout) :: x2
      real*8 , optional, dimension(*), intent(inout) :: x3
      real*8 , optional, dimension(*), intent(inout) :: x4
      real*8 , optional,               intent(inout) :: energy
      integer, optional,               intent(in)    :: ndim
      integer, optional,               intent(in)    :: mdim
      integer, optional,               intent(in)    :: odim
      integer, optional,               intent(in)    :: pdim
      integer, optional,               intent(in)    :: state
      integer, optional,               intent(in)    :: stateL
      integer, optional,               intent(in)    :: msproj
      integer, optional,               intent(in)    :: msprojL
      integer, optional,               intent(in)    :: multiplet
      integer, optional,               intent(in)    :: multipletL
      logical, optional,               intent(in)    :: rdm1
      logical, optional,               intent(in)    :: rdm2
      logical, optional,               intent(in)    :: rdm3
      logical, optional,               intent(in)    :: rdm4
      logical, optional,               intent(in)    :: Key_CION
      integer, optional,               intent(in)    :: iterSCF
      ! Leon 02-12-2016: added optional custom checkpoint names
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2

      !print *, 'DMRG interface called with task == ',trim(task)

      if(dmrg_host_program_settings%myrank == 0)then

        call qcmaquis_interface_ctl(                                  &
                                    task,                             &
                                    x1,                               &
                                    x2,                               &
                                    x3,                               &
                                    x4,                               &
                                    energy,                           &
                                    ndim,                             &
                                    mdim,                             &
                                    odim,                             &
                                    pdim,                             &
                                    state,                            &
                                    stateL,                           &
                                    msproj,                           &
                                    msprojL,                          &
                                    multiplet,                        &
                                    multipletL,                       &
                                    rdm1,                             &
                                    rdm2,                             &
                                    rdm3,                             &
                                    rdm4,                             &
                                    Key_CION,                         &
                                    IterSCF,                          &
                                    checkpoint1,                      &
                                    checkpoint2                       &
                                   )

      end if ! only myrank == 0 enters the interface

      !> update all processes if needed
      if(dmrg_host_program_settings%nprocs > 1 .and. dmrg_host_program_settings%runs_parallel)then
        call dmrg_task_process_update(                 &
                                      task   = task,   &
                                      x1     = x1,     &
                                      x2     = x2,     &
                                      energy = energy, &
                                      ndim   = ndim,   &
                                      mdim   = mdim,   &
                                      rdm1   = rdm1,   &
                                      rdm2   = rdm2    &
                                     )
      end if

  end subroutine dmrg_interface_ctl

  subroutine dmrg_task_process_update(                                  &
                                      task,                             &
                                      x1,                               &
                                      x2,                               &
                                      energy,                           &
                                      ndim,                             &
                                      mdim,                             &
                                      rdm1,                             &
                                      rdm2                              &
                                     )
      !> purpose: either update all co-worker processes or let them wait until the master finishes a DMRG-related task

      character(len=8),                intent(in)    :: task
      real*8 , optional, dimension(*), intent(inout) :: x1
      real*8 , optional, dimension(*), intent(inout) :: x2
      real*8 , optional,               intent(inout) :: energy
      integer, optional,               intent(in)    :: ndim
      integer, optional,               intent(in)    :: mdim
      logical, optional,               intent(in)    :: rdm1
      logical, optional,               intent(in)    :: rdm2
      integer                                        :: error

#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"

      ! select task
      select case(trim(task))
        case('fci dump')
           call GA_sync()
        case('overlap ', 'overlapR')
           Call GA_Brdcst(MT_DBL, energy, 1, 0)
        case('imp spdX', 'imp rdmY')
           Call GA_Brdcst(MT_DBL, x1, nDim, 0)
        case('imp rdmX')
           !if(present(rdm1) .and. rdm1) print *, 'rdm1 before ',x1(1:10)
           if(present(rdm1) .and. rdm1) Call GA_Brdcst(MT_DBL, x1, nDim, 0)
           if(present(rdm2) .and. rdm2) Call GA_Brdcst(MT_DBL, x2, mDim, 0)
           call GA_sync()
           !if(present(rdm1) .and. rdm1) print *, 'rdm1 after ',x1(1:10)
        case('run DMRG')
           call GA_sync()
           Call GA_Brdcst(MT_DBL, dmrg_energy%dmrg, 1, 0)
           Call GA_Brdcst(MT_DBL, dmrg_energy%dmrg_state_specific, size(dmrg_energy%dmrg_state_specific), 0)
           call GA_sync()
        case default
          write(6,*) 'nothing to synchronize...'
          call Abend()
      end select

#endif

      ! take care of unused variable warnings
      if (.false.) then
        call Unused_integer(ndim)
        call Unused_integer(mdim)
        call Unused_integer(error)
        call Unused_logical(rdm1)
        call Unused_logical(rdm2)
        call Unused_character(task)
        call Unused_real(energy)
        call Unused_real_array(x1)
        call Unused_real_array(x2)
      end if

      end subroutine dmrg_task_process_update

#else
implicit none
integer :: qcmaquis_wrapper_dummy
#endif

end module qcmaquis_interface_wrapper
