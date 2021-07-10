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
subroutine prepare_decay
  use rhodyn_data
  use rhodyn_utils, only: transform, mult, dashes
  use definitions, only: u6
  implicit none
! energy and time should fullfill relation delta_E*delta_t=h
! construct the decay in SOC states basis sets

  if (ipglob>3) write(u6,*) 'Begin of prepare_decay'

  decay=zero
! Auger decay
  if (flag_decay) then
    do i=Nval+1,Nval+N_L3
      decay(i,i)=-tau_L3/2/pi
    enddo
    do i=Nval+N_L3+1,Nstate
      decay(i,i)=-tau_L2/2/pi
    enddo
    if (basis=='CSF') then
      call mult(CSF2SO,decay,tmp)
      call mult(tmp,CSF2SO,decay,.False.,.True.)
    elseif (basis=='SF') then
      call mult(SO_CI,decay,tmp)
      call mult(tmp,SO_CI,decay,.False.,.True.)
    endif
  endif

! ionization
  if (flag_dyson.and.ion_diss/=0d0) then
    ii=1
    do k=1,N
      do i=ii,(ii+nconf(k)*ispin(k)-1)
        if (ion_blocks(k)) decay(i,i) = decay(i,i) - ion_diss
      enddo
      ii=ii+nconf(k)*ispin(k)
    enddo
    if (basis=='CSF') then
      call mult(U_CI_compl,decay,tmp)
      call mult(tmp,U_CI_compl,decay,.False.,.True.)
    elseif (basis=='SO') then
      call mult(SO_CI,decay,tmp,.True.,.False.)
      call mult(tmp,SO_CI,decay)
    endif
  endif

!!!!!!!!!!
  if (ipglob>4) then
    call dashes()
    write(u6,*) 'Decay matrix'
    do i=1,Nstate
      write(u6,*)(decay(i,j),j=1,Nstate)
    enddo
    call dashes()
  endif
!!!!!!!!!!

  if (ipglob>3) write(u6,*) 'End of prepare_decay'

end
