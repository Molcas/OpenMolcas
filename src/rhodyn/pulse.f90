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
  use rhodyn_data, only: zero, i, E_field, N_pulse, pulse_type,&
                         pulse_vec, pulse_vector, amp, taushift,&
                         sigma, power_shape, omega, phi,&
                         lu_pls, temp_vec, out_t, out_pulse,&
                         dipole_basis
  use mh5, only: mh5_put_dset
  use constants, only: pi, auToFs
  implicit none
  complex(8), dimension(:,:), intent(in) :: H0
  complex(8), dimension(:,:), intent(out) :: Ht
  real(8), intent(in)       :: time
  integer, intent(in), optional :: count

  E_field=zero

  do i=1,N_pulse
!
! sine^N pulse
! add description here
    if (pulse_type(1:3)=='SIN') then
      pulse_vec(:)=pulse_vector(i,:)
      if (abs(time-taushift(i))<=sigma(i)) then
        E_field = amp(i)*pulse_vec &
        *sin( pi*(time-taushift(i)) / (2.d0*sigma(i)) )**power_shape &
        *sin(omega(i) * (time-taushift(i)) + phi(i))
! add correction here
      else
        E_field=0d0
      endif
!
! cosine^N pulse
! add description here
    elseif (pulse_type(1:3)=='COS') then
      pulse_vec(:)=pulse_vector(i,:)
      if (abs(time-taushift(i))<=sigma(i)) then
        E_field = amp(i)*pulse_vec &
        *cos( pi*(time-taushift(i)) / (2.d0*sigma(i)) )**power_shape &
        *sin(omega(i) * (time-taushift(i)) + phi(i))
! add correction here
      else
        E_field=0d0
      endif
!
! gaussian pulse
! add description here
    elseif (pulse_type=='GAUSS') then
      pulse_vec(:)=pulse_vector(i,:)
      E_field = amp(i)*pulse_vec &
      *exp(-(time-taushift(i))**2 / (2*sigma(i)**2)) &
      *sin(omega(i) * (time-taushift(i)) + phi(i))
      ! add correction here
!
! monochromatic pulse
! add description here
    elseif (pulse_type=='MONO') then
      pulse_vec(:)=pulse_vector(i,:)
      E_field = amp(i)*pulse_vec &
      *sin(omega(i) * (time-taushift(i)) + phi(i))
!
! train of gaussian pulses
! elseif (pulse_type=='Train_Diff') then
!   E_field = E_field+amp(i)*pulse_vector(i,:)* &
!             exp(-(time-(tau+(i-1)*shift(i)))**2/ &
!             (2*sigma(i)**2))*sin(omega(i)*time+phi(i))
!
! explicitely polarized pulses
! think of more clever definition
    elseif (pulse_type=='MONO_R_CIRCLE') then
      E_field(1) = amp(i)*sin(omega(i)*time+phi(i))
      E_field(2) = Amp(1)*cos(omega(i)*time+phi(i))
      E_field(3) = 0d0
!
    elseif (pulse_type=='MONO_L_CIRCLE') then
      E_field(1) = amp(i)*cos(omega(i)*time+phi(i))
      E_field(2) = amp(i)*sin(omega(i)*time+phi(i))
      E_field(3) = 0d0
!
    elseif (pulse_type=='GAUSS_R_CIRCLE') then
      E_field(1) = amp(i)*sin(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(2) = amp(i)*cos(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(3)=0d0
!
    elseif (pulse_type=='GAUSS_L_CIRCLE')then
      E_field(1) = amp(i)*cos(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(2) = amp(i)*sin(omega(i)*time+phi(i))* &
                   exp(-(time-taushift(i))**2/(2*sigma(i)**2))
      E_field(3) = 0d0
    endif
  enddo
!
! saving pulse to files
  if (present(count)) then
    write(lu_pls,'(7(G25.15E3,2X))') time*auToFs, &
                  (dble(E_field(i)),aimag(E_field(i)), i=1,3)
    do i=1,3
      temp_vec(2*i-1) = dble (E_field(i))
      temp_vec(2*i)   = aimag(E_field(i))
    enddo
    call mh5_put_dset(out_pulse,temp_vec, [1,6],[count-1,0])
    call mh5_put_dset(out_t, [time*auToFs], [1], [count-1])
  endif
!
! update Hamiltonian H0 to Ht adding
! scalar product of E_field and dipole moment
  Ht = H0
  do i=1,3
    Ht(:,:) = Ht(:,:) + dipole_basis(:,:,i)*E_field(i)
  enddo
!
end
