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
* Copyright (C) 2020, Giovanni Li Manni                                *
************************************************************************
      module print_RDMs_NECI_format
      use index_symmetry, only : two_el_idx, one_el_idx
      implicit none
      private
      public:: printRDMs_NECI
      contains

      subroutine printRDMs_NECI(DMAT,NAC,PMAT,PA,NACPAR)
      integer, intent(in) :: NAC, NACPAR
      real*8, intent(in)  :: DMAT(NAC), PMAT(NACPAR), PA(NACPAR)
      real*8, parameter :: thrsh = 1.0d-12
      integer :: i, t_idx, u_idx, v_idx, x_idx

        Write(6,*) ' In printRDMs_NECI:'

        do i = 1, (NACPAR*(NACPAR+1)/2)
          call two_el_idx(i, t_idx, u_idx, v_idx, x_idx)
          if(v_idx /= x_idx) then
             if(abs(PMAT(i)+PA(i)).gt.thrsh) then
                 write(6,'(1X,4I5,F20.12)')
     &           t_idx, u_idx, v_idx, x_idx, (PMAT(i)+PA(i))
             end if
             if (abs(PMAT(i)-PA(i)).gt.thrsh) then
                 write(6,'(1X,4I5,F20.12)')
     &           t_idx, u_idx, x_idx, v_idx, (PMAT(i)-PA(i))
             end if
          else
             if (abs(PMAT(i)*2.0d0).gt.thrsh) then
                write(6,'(1X,4I5,F20.12)')
     &          t_idx, u_idx, v_idx, x_idx, PMAT(i)*2.0d0
             end if
          end if
        end do

        do i = 1, (NAC*(NAC+1)/2)
          call one_el_idx(i, t_idx, u_idx)
          if(abs(DMAT(i)).gt.thrsh) then
             write(6,'(1X,4I5,F20.12)')
     &          t_idx, u_idx, 0, 0, DMAT(i)
          end if
        end do

        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',DMAT,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',PMAT,NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',PA , NACPAR)

      end subroutine printRDMs_NECI

      end module print_RDMs_NECI_format
