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

subroutine prepare_decay()
! energy and time should fullfill relation delta_E*delta_t=h
! construct the decay in SOC states basis sets

use rhodyn_data, only: basis, CSF2SO, decay, flag_decay, flag_dyson, ion_diss, ipglob, ispin, N, N_L3, nconf, Nstate, Nval, SO_CI, &
                       tau_L2, tau_L3, tmp, U_CI_compl
use rhodyn_utils, only: mult, dashes
use Constants, only: Zero, cZero, pi
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, j, k, ii
logical(kind=iwp), parameter :: ion_blocks(5) = [.true.,.false.,.true.,.false.,.true.]

if (ipglob > 3) write(u6,*) 'Begin of prepare_decay'

decay(:,:) = cZero
! Auger decay
if (flag_decay) then
  do i=Nval+1,Nval+N_L3
    decay(i,i) = -tau_L3/2/pi
  end do
  do i=Nval+N_L3+1,Nstate
    decay(i,i) = -tau_L2/2/pi
  end do
  if (basis == 'CSF') then
    call mult(CSF2SO,decay,tmp)
    call mult(tmp,CSF2SO,decay,.false.,.true.)
  else if (basis == 'SF') then
    call mult(SO_CI,decay,tmp)
    call mult(tmp,SO_CI,decay,.false.,.true.)
  end if
end if

! ionization
if (flag_dyson .and. (ion_diss /= Zero)) then
  ii = 1
  do k=1,N
    do i=ii,(ii+nconf(k)*ispin(k)-1)
      if (ion_blocks(k)) decay(i,i) = decay(i,i)-ion_diss
    end do
    ii = ii+nconf(k)*ispin(k)
  end do
  if (basis == 'CSF') then
    call mult(U_CI_compl,decay,tmp)
    call mult(tmp,U_CI_compl,decay,.false.,.true.)
  else if (basis == 'SO') then
    call mult(SO_CI,decay,tmp,.true.,.false.)
    call mult(tmp,SO_CI,decay)
  end if
end if

!!!!!!!!!!
if (ipglob > 4) then
  call dashes()
  write(u6,*) 'Decay matrix'
  do i=1,Nstate
    write(u6,*) (decay(i,j),j=1,Nstate)
  end do
  call dashes()
end if
!!!!!!!!!!

if (ipglob > 3) write(u6,*) 'End of prepare_decay'

end subroutine prepare_decay
