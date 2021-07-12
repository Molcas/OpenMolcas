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
subroutine rhodyn_init()
!***********************************************************************
!
!     input parameters set by default. Later they can be overwritten
!     when reading input file
!
!***********************************************************************
  use rhodyn_data
  use stdalloc, only: mma_allocate, mma_deallocate
  use constants, only: auToFs, auToCm, auToeV
  implicit none

! preparation 1 means standard program workflow (see rhodyn_data)
  preparation   = 1
  flag_test     = .False.
  p_style       = 'SF'
  N_Populated   = 1
! temperature T needed if p_style includes 'thermal'
  T             = 300
  Nmode         = 0
! be default propagation basis is spin free states
  basis         = 'SF'
  tout          = 0.05d0/auToFs
  initialtime   = 0.0d0/auToFs
  finaltime     = 10.0d0/auToFs
  timestep      = 0.0005d0/auToFs
  Method        = 'classic_RK4'
  errorthreshold= 1.0d-06
  safety        = 0.9
  deltaE        = 50d0/auToCm
  V             = 100d0/auToCm
  Nval          = 160
  N_L3          = 175
  tau_L3        = 0.4d0/autoev
  N_L2          = 585
  tau_L2        = 1.04d0/autoev
  flag_dyson    = .False.
  alpha         = 1d-3
  ion_diss      = 0d0
  ion_blocks    = (/.True.,.False.,.True.,.False.,.True./)
  flag_diss     = .False.
  gamma         = 300/auToCm
  HRSO          = .False.
  kext          = .False.
  DM_basis      = 'SF_SO'
! full density matrix saving time step
  time_fdm      = 1.0d0/auToFs
! general idea is that additional features are disabled by default
! except for pulse flag
  flag_so       = .False.
  flag_decay    = .False.
  flag_fdm      = .False.
  flag_dipole   = .False.
  flag_emiss    = .False.
  flag_pulse    = .True.
  Pulse_type    = 'Gaussian'
! number of incoming pulses N_pulse supposed to be 1 by default
! later when reading input it can be changed with
! reallocation of all corresponding arrays
  N_pulse       = 1
  call mma_allocate(amp,N_pulse)
  call mma_allocate(taushift,N_pulse)
  call mma_allocate(pulse_vector,N_pulse,3)
  call mma_allocate(sigma,N_pulse)
  call mma_allocate(omega,N_pulse)
  call mma_allocate(phi,N_pulse)
  amp(1)           = 2.5d0
  taushift(1)      = 3.0d0/auToFs
  pulse_vector(1,1)= one
  pulse_vector(1,2)= zero
  pulse_vector(1,3)= zero
  sigma(1)         = 1.0d0/auToFs
  omega(1)         = 710d0/autoev
  phi(1)           = 0d0

end
