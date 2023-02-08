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
! Copyright (C) 2021,2023, Vladislav Kochetov                          *
!***********************************************************************

subroutine pulse(H0,Ht,time,pcount)
!***********************************************************************
! Purpose :  construct the time-dependent hamiltonian
!                        Ht = H0 + dipole*E_field(time)
!            argument 'pcount' is for storage of pulse
!            if pcount=-1 then current value of pulse is not stored
!***********************************************************************

use rhodyn_data, only: amp, d, dipole_basis, flag_acorrection, linear_chirp, lu_pls, N_pulse, omega, out_pulse, out_t, phi, &
                       power_shape, pulse_type, pulse_vec, pulse_vector, sigma, taushift
use mh5, only: mh5_put_dset
use Constants, only: Two, cZero, pi, auToFs
use Definitions, only: wp, iwp

implicit none
complex(kind=wp), intent(in) :: H0(d,d)
complex(kind=wp), intent(out) :: Ht(d,d)
real(kind=wp), intent(in) :: time
integer(kind=iwp), intent(in) :: pcount
integer(kind=iwp) :: i
real(kind=wp) :: omega_local, t_local, temp_vec(6)
complex(kind=wp) :: E_field(3)

E_field(:) = cZero

do i=1,N_pulse
  t_local = time-taushift(i)
  omega_local = omega(i)+linear_chirp*t_local
  pulse_vec(:) = pulse_vector(i,:)

  ! sine^N pulse
  ! A\vec{e}\sin^n(\pi(t-t_0)/(2\sigma))\sin{(\Omega(t-t_0)+\varphi_0)}
  ! duration of the pulse equals 2\sigma
  if (pulse_type(1:3) == 'SIN') then
    if (abs(t_local) <= sigma(i)) then
      E_field = E_field+amp(i)*pulse_vec*sin(pi*t_local/(Two*sigma(i)))**power_shape*sin(omega_local*t_local+phi(i))
      ! correction due to vector potential derivative [Paramonov_JPCA_2012]:
      ! -A\vec{e}\pi/(2\sigma\Omega)
      !     \sin(\pi(t-t_0)/(\sigma))\cos{(\Omega(t-t_0)+\varphi_0)}
      ! only for sine square shape form
      if (flag_acorrection .and. (power_shape == 2)) then
        E_field = E_field-amp(i)*pulse_vec*pi/(2*sigma(i)*omega_local)*sin(pi*t_local/sigma(i))*cos(omega_local*t_local+phi(i))
      !else
      !  E_field = Zero
      end if
    end if

  ! cos^N pulse:
  ! A\vec{e}\cos^n(\pi(t-t_0)/(2\sigma))\sin{(\Omega(t-t_0)+\varphi_0)}
  ! duration of the pulse equals 2\sigma
  else if (pulse_type(1:3) == 'COS') then
    if (abs(t_local) <= sigma(i)) then
      E_field = E_field+amp(i)*pulse_vec*cos(pi*t_local/(Two*sigma(i)))**power_shape*sin(omega_local*t_local+phi(i))
      ! correction is the same as in sine case, but with different choice of
      ! vector potential, here vector potential is given:
      ! -\vec{e}\frac{A}{\Omega} \cos^2(\pi(t-t_0)/(2\sigma))
      !     \cos{(\Omega(t-t_0)+\varphi_0)}
      ! only for cosine square shape form
      if (flag_acorrection .and. (power_shape == 2)) then
        E_field = E_field-amp(i)*pulse_vec*pi/(2*sigma(i)*omega_local)*sin(pi*t_local/sigma(i))*cos(omega_local*t_local+phi(i))
      !else
      !  E_field = Zero
      end if
    end if

  ! gaussian pulse
  ! A\vec{e}exp{-(t-t_0)^2/(2\sigma^2)}\sin{(\Omega(t-t_0)+\varphi_0)}
  else if (pulse_type == 'GAUSS') then
    E_field = E_field+amp(i)*pulse_vec*exp(-t_local**2/(2*sigma(i)**2))*sin(omega_local*t_local+phi(i))
    ! correction due to vector potential derivative:
    ! \vec{e}\frac{A(t-t_0)}{\sigma^2\Omega}
    !     \exp{-(t-t_0)^2/(2\sigma)^2}\cos{(\Omega(t-t_0) + \varphi_0)}
    if (flag_acorrection) then
      E_field = E_field+amp(i)*pulse_vec*t_local/(sigma(i)**2*omega_local)*exp(-t_local**2/(2*sigma(i)**2))* &
                cos(omega_local*t_local+phi(i))
    end if

  ! monochromatic pulse
  ! A\vec{e}\sin{(\Omega(t-t_0)+\varphi_0)}
  else if (pulse_type == 'MONO') then
    E_field = amp(i)*pulse_vec*sin(omega_local*(time-taushift(i))+phi(i))

  ! explicitly polarized pulses
  ! think of more clever definition
  else if (pulse_type == 'MONO_R_CIRCLE') then
    E_field(1) = amp(i)*sin(omega(i)*time+phi(i))
    E_field(2) = Amp(1)*cos(omega(i)*time+phi(i))
    E_field(3) = cZero

  else if (pulse_type == 'MONO_L_CIRCLE') then
    E_field(1) = amp(i)*cos(omega(i)*time+phi(i))
    E_field(2) = amp(i)*sin(omega(i)*time+phi(i))
    E_field(3) = cZero

  else if (pulse_type == 'GAUSS_R_CIRCLE') then
    E_field(1) = amp(i)*sin(omega(i)*time+phi(i))*exp(-(time-taushift(i))**2/(Two*sigma(i)**2))
    E_field(2) = amp(i)*cos(omega(i)*time+phi(i))*exp(-(time-taushift(i))**2/(Two*sigma(i)**2))
    E_field(3) = cZero

  else if (pulse_type == 'GAUSS_L_CIRCLE') then
    E_field(1) = amp(i)*cos(omega(i)*time+phi(i))*exp(-(time-taushift(i))**2/(Two*sigma(i)**2))
    E_field(2) = amp(i)*sin(omega(i)*time+phi(i))*exp(-(time-taushift(i))**2/(Two*sigma(i)**2))
    E_field(3) = cZero
  end if
end do

! saving pulse to files
if (pcount >= 0) then
  write(lu_pls,'(7(g25.15e3,2x))') time*auToFs,(real(E_field(i)),aimag(E_field(i)),i=1,3)
  do i=1,3
    temp_vec(2*i-1) = real(E_field(i))
    temp_vec(2*i) = aimag(E_field(i))
  end do
  call mh5_put_dset(out_pulse,temp_vec,[1,6],[pcount-1,0])
  call mh5_put_dset(out_t,[time*auToFs],[1],[pcount-1])
end if

! update Hamiltonian H0 to Ht adding
! scalar product of E_field and dipole moment
Ht = H0
do i=1,3
  Ht(:,:) = Ht(:,:)+dipole_basis(:,:,i)*E_field(i)
end do

end subroutine pulse
