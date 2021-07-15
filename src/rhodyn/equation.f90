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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine equation (time,rhot,res)
  use rhodyn_data, only: pulse_func, flag_pulse, d, onei, zero, one, &
                         flag_decay, ion_diss, flag_diss, decay, &
                         hamiltonian, hamiltoniant, &
                         K_bar_basis, kab_basis, i, j
  implicit none
!
!***********************************************************************
! Purpose : Liouville equation is solved here
!
!***********************************************************************
!
!  time   : current time
!  rhot   : density matrix at current time
!  res    : obtained left part of Liouville equation d(rhot)/d(time)
!
  real(8), intent(in) :: time
  complex(8), dimension(:,:), intent(in) :: rhot
  complex(8), dimension(:,:), intent(out):: res
  procedure(pulse_func) :: pulse

! if pulse in enableb, modify hamiltonian at time t:
  if (flag_pulse) call pulse(hamiltonian,hamiltoniant,time)

! get right part of Liouville equation
  call zgemm_('N','N',d,d,d,-onei,hamiltoniant,d,rhot,d,zero,res,d)
  call zgemm_('N','N',d,d,d, onei,rhot,d,hamiltoniant,d, one,res,d)

! auger decay part
  if (flag_decay.or.ion_diss/=0d0) then
    call zgemm_('N','N',d,d,d,one,decay,d,rhot,d,one,res,d)
  endif

! if dissipation (nuclear bath) is considered
  if (flag_diss) then
    do i=1,d
      do j=1,d
        if (i/=j) then
          res(i,j) = res(i,j) - K_bar_basis(i,j) * rhot(i,j)
        endif
        res(i,i) = res(i,i) - Kab_basis(i,j) * rhot(i,i) + &
                              Kab_basis(j,i) * rhot(j,j)
      enddo
    enddo
  endif

end
