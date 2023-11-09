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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine stat2_cvb()

use casvb_global, only: cpu0, cpu_prev, endvar, ipr, n_2el, n_applyh, n_applyt, n_cihess, n_hess, n_orbhess
use Definitions, only: wp, u6

implicit none
real(kind=wp), external :: tim_cvb

if (ipr(3) >= 1) then
  write(u6,'(/,a,i16)') ' Total number of structure transformations :',n_applyt
  write(u6,'(a,i17)') ' Total number of Hamiltonian applications :',n_applyh
  write(u6,'(a,i11)') ' Total number of 2-electron density evaluations :',n_2el
  write(u6,'(a,i21)') ' Total number of Hessian applications :',n_hess
  if (n_orbhess > 0) write(u6,'(a,i8)') ' Total number of pure orbital Hessian applications :',n_orbhess
  if (n_cihess > 0) write(u6,'(a,i13)') ' Total number of pure CI Hessian applications :',n_cihess
  if (.not. endvar) then
    call date2_cvb(tim_cvb(cpu0))
  else
    call date2_cvb(tim_cvb(cpu0)+cpu_prev)
  end if
else
  cpu_prev = cpu_prev+tim_cvb(cpu0)
end if

return

end subroutine stat2_cvb
