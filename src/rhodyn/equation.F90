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

subroutine equation(time,rhot,res)
!***********************************************************************
! Purpose : RHS of Liouville equation is obtained here
!
!***********************************************************************
!
!  time   : current time
!  rhot   : density matrix at current time
!  res    : obtained RHS of Liouville equation d(rhot)/d(time)

use rhodyn_data, only: d, decay, flag_decay, flag_diss, flag_pulse, hamiltonian, hamiltoniant, ion_diss, K_bar_basis, kab_basis
!use rhodyn_utils, only: print_c_matrix
use Constants, only: Zero, cZero, cOne, Onei
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: time
complex(kind=wp), intent(in) :: rhot(d,d)
complex(kind=wp), intent(out) :: res(d,d)
integer(kind=iwp) :: i, j
!procedure(pulse_func) :: pulse

! if pulse is enabled, modify Hamiltonian at time t:
if (flag_pulse) call pulse(hamiltonian,hamiltoniant,time,-1)

!!! debug !!!
!call print_c_matrix(hamiltoniant, d, 'hamiltoniant', 6)
!call print_c_matrix(rhot, d, 'rhot', 6)

! get right part of Liouville equation -i*(hamiltoniant*rhot - rhot*hamiltoniant)
call zgemm_('N','N',d,d,d,-Onei,hamiltoniant,d,rhot,d,cZero,res,d)
call zgemm_('N','N',d,d,d,Onei,rhot,d,hamiltoniant,d,cOne,res,d)

!call print_c_matrix(res, d, 'res', 6)

! Auger decay part
if (flag_decay .or. (ion_diss /= Zero)) then
  call zgemm_('N','N',d,d,d,cOne,decay,d,rhot,d,cOne,res,d)
end if

! if dissipation (nuclear bath) is considered
if (flag_diss) then
  do i=1,d
    do j=1,d
      if (i /= j) then
        res(i,j) = res(i,j)-K_bar_basis(i,j)*rhot(i,j)
      end if
      res(i,i) = res(i,i)-Kab_basis(i,j)*rhot(i,i)+Kab_basis(j,i)*rhot(j,j)
    end do
  end do
end if

end subroutine equation
