************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Oskar Weser                                      *
************************************************************************
      module CI_solver_util
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      implicit none
      private
      public :: wait_and_read, abort_
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
      integer*4 :: error
      integer*4, parameter :: one4=1, root4=0
#endif

      interface
        integer function isfreeunit(iseed)
          integer, intent(in) :: iseed
        end function
      end interface

      contains


      subroutine wait_and_read(energy)
        real*8, intent(out) :: energy
        logical :: newcycle_found
        integer :: LuNewC
        newcycle_found = .false.
        do while(.not. newcycle_found)
          call sleep(1)
          if (myrank == 0) call f_Inquire('NEWCYCLE', newcycle_found)
#ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call MPI_Bcast(newcycle_found, one4, MPI_LOGICAL,
     &                     root4, MPI_COMM_WORLD, error)
          end if
#endif
        end do
        if (myrank == 0) then
          write(6, *) 'NEWCYCLE file found. Proceding with SuperCI'
          LuNewC = isFreeUnit(12)
          call molcas_open(LuNewC, 'NEWCYCLE')
            read(LuNewC,*) energy
          close(LuNewC, status='delete')
          write(6, *) 'I read the following energy:', energy
        end if
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          call MPI_Bcast(energy, one4, MPI_REAL8,
     &                   root4, MPI_COMM_WORLD, error)
        end if
#endif
      end subroutine wait_and_read

      subroutine abort_(message)
        character(*), intent(in) :: message
        call WarningMessage(2, message)
        call QTrace()
        call Abend()
      end subroutine



      end module CI_solver_util
