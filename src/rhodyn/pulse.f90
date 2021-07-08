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
subroutine pulse(H0,Ht,time,count)
!***********************************************************************
!
! Purpose :  construct the time-dependent hamiltonian
!                                       H_field=dipole*E_field
!
!***********************************************************************
  use rhodyn_data
  use mh5
  implicit none
  complex(8),dimension(:,:) :: H0,Ht
  real(8)                   :: time
  integer,optional          :: count

  Ht = zero
  E_field=zero

  if (pulse_type=='Sin2') then
    pulse_vec(:)=pulse_vector(1,:)
! the time scale only vaild for sigma=5.0
    if (abs(time-tau)<=sigma(1)) then
      E_field = Amp(1)*pulse_vec*(sin(pi*time/ &
              Sigma(1)/2))**2*sin(omega(1)*time+phi(1))
    else
      E_field=0d0
    endif

  elseif (pulse_type=='Cos2') then
    pulse_vec(:)=pulse_vector(1,:)
    if (abs(time-tau)<=sigma(1)*pi/2.d0) then
      E_field = amp(1)*pulse_vec &
                *cos((time-tau)/sigma(1))**2 &
                *cos(omega(1)*(time-tau)+phi(1))
    else
      E_field=0d0
    endif

  elseif (pulse_type=='Gaussian') then
    pulse_vec(:)=pulse_vector(1,:)
    E_field = amp(1)*pulse_vec &
              *exp(-(time-tau)**2/(2*sigma(1)**2)) &
              *sin(omega(1)*(time-tau)+phi(1))

  elseif (Pulse_type=='Continue') then
    pulse_vec(:)=Pulse_vector(1,:)
    E_field=Amp(1)*pulse_vec*sin(omega(1)*time+phi(1))

  elseif (pulse_type=='Train_Same') then
! in this case the lights are same,
! so we chose the data in the first dimension
    pulse_vec(:)=pulse_vector(1,:)
    do i=1,N_pulse
      E_field = E_field+amp(1)*pulse_vec*exp(-(time-(tau+i* &
                shift(2)-shift(2)))**2/(2*Sigma(1)**2))* &
                sin(omega(1)*time+phi(1))
    enddo

  elseif (pulse_type=='Train_Diff') then
    do i=1,N_pulse
      E_field = E_field+amp(i)*pulse_vector(i,:)* &
                exp(-(time-(tau+(i-1)*shift(i)))**2/ &
                (2*sigma(i)**2))*sin(omega(i)*time+phi(i))
    enddo

  elseif (pulse_type=='R_Circle_Continue') then
        E_field(1) = Amp(1)*sin(omega(1)*time+phi(1))
        E_field(2) = Amp(1)*cos(omega(1)*time+phi(1))
        E_field(3) = 0d0

  elseif (pulse_type=='L_Circle_Continue') then
        E_field(1) = Amp(1)*cos(omega(1)*time+phi(1))
        E_field(2) = Amp(1)*sin(omega(1)*time+phi(1))
        E_field(3) = 0d0

  elseif (pulse_type=='R_Circle_Gaussian') then
        E_field(1) = Amp(1)*sin(omega(1)*time+phi(1))* &
                 exp(-(time-tau)**2/(2*Sigma(1)**2))
        E_field(2) = Amp(1)*cos(omega(1)*time+phi(1))* &
                 exp(-(time-tau)**2/(2*Sigma(1)**2))
        E_field(3)=0d0

  elseif (pulse_type=='L_Circle_Gaussian')then
        E_field(1) = Amp(1)*cos(omega(1)*time+phi(1))* &
                 exp(-(time-tau)**2/(2*Sigma(1)**2))
        E_field(2) = Amp(1)*sin(omega(1)*time+phi(1))* &
                 exp(-(time-tau)**2/(2*Sigma(1)**2))
        E_field(3) = 0d0

  elseif (pulse_type=='newL_Circle_Continue') then
        E_field(1)=Amp(1)*sin(omega(1)*time+phi(1))
        E_field(2)=-Amp(1)*cos(omega(1)*time+phi(1))
        E_field(3)=0d0
  endif

  if (present(count)) then
    write(lu_pls,'(7(G25.15E3,2X))') time/fstoau, &
                  (dble(E_field(i)),aimag(E_field(i)), i=1,3)
    do i=1,3
      temp_vec(2*i-1) = dble (E_field(i))
      temp_vec(2*i)   = aimag(E_field(i))
    enddo
    call mh5_put_dset(out_pulse,temp_vec, [1,6],[count-1,0])
    call mh5_put_dset(out_t, [time/fstoau], [1], [count-1])
  endif

  do i=1,3
    Ht(:,:)=Ht(:,:)+dipole_basis(:,:,i)*E_field(i)
  enddo
  Ht = H0 + Ht

end
