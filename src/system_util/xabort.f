************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine xabort(rc)
C     this routine aborts the process(es) with rc
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      implicit none
      integer :: rc
#ifdef _MOLCAS_MPP_
#  include "para_info.fh"
      integer*4 :: rc4, ierr4
#endif

#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
C     unline abort, mpi_abort does not print a backtrace,
C     so we do it manually here
        call xbacktrace
        rc4 = int(rc,kind(rc4))
        call mpi_abort(MPI_COMM_WORLD,rc4,ierr4)
      else
#endif

#if defined (__GNUC__) || defined(__INTEL_COMPILER)
        call abort
#else
        stop 'Molcas aborted...'
#endif

#ifdef _MOLCAS_MPP_
      end if
#else
c Avoid unused argument warnings
      if (.false.) call Unused_integer(rc)
#endif
      end
