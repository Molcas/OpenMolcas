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

subroutine xabort(rc)
! this routine aborts the process(es) with rc

#ifdef _MOLCAS_MPP_
use mpi, only: MPI_COMM_WORLD
use Para_Info, only: Is_Real_Par
use Definitions, only: MPIInt
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: rc
#ifdef _MOLCAS_MPP_
integer(kind=MPIInt) :: rc4, ierr4
#else
#include "macros.fh"
unused_var(rc)
#endif

#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  ! unlike abort, mpi_abort does not print a backtrace,
  ! so we do it manually here
  call xbacktrace()
  rc4 = int(rc,kind=MPIInt)
  call mpi_abort(MPI_COMM_WORLD,rc4,ierr4)
else
#endif

#if defined (__GNUC__) || defined(__INTEL_COMPILER)
  call abort()
#else
  stop 'Molcas aborted...'
#endif

#ifdef _MOLCAS_MPP_
end if
#endif

return

end subroutine xabort
