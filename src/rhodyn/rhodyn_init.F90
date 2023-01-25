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
! input parameters set by default. Later they can be overwritten
! when reading input file
!***********************************************************************

use Constants, only: Zero, One, Three, Ten, cZero, cOne, auToFs, auToCm, auToeV
use Definitions, only: wp
use rhodyn_data, only: alpha, amp, basis, DM_basis, errorthreshold, finaltime, flag_acorrection, flag_decay, flag_dipole, &
                       flag_diss, flag_dyson, flag_emiss, flag_fdm, flag_pulse, flag_so, flag_test, cgamma, HRSO, initialtime, &
                       ion_diss, kext, k_max, linear_chirp, method, N_L2, N_L3, N_Populated, N_pulse, Nmode, Nval, omega, p_style, &
                       phi, pulse_type, pulse_vector, runmode, safety, sigma, T, tau_L2, tau_L3, taushift, time_fdm, timestep, tout
use stdalloc, only: mma_allocate

implicit none

! runmode 1 means standard program workflow (see rhodyn_data)
runmode = 1
flag_test = .false.
p_style = 'SF'
N_Populated = 1
! temperature T needed if p_style includes 'thermal'
T = 300
Nmode = 0
! be default propagation basis is spin free states
basis = 'SF'
DM_basis = 'SF_SO'
tout = 0.05_wp/auToFs
initialtime = Zero/auToFs
finaltime = Ten/auToFs
timestep = 0.0005_wp/auToFs
method = 'CLASSIC_RK4'
errorthreshold = 1.0e-06_wp
! safety parameter for adaptive-size methods can be set to 0.95
! for acceleration of calculations
safety = 0.9_wp
!deltaE = 50.0_wp/auToCm
!V = 100.0_wp/auToCm
k_max = 2
Nval = 160
N_L3 = 175
tau_L3 = 0.4_wp/autoev ! Auger decay rate for Fe L3
N_L2 = 585
tau_L2 = 1.04_wp/autoev ! Auger decay rate for Fe L2
flag_dyson = .false.
alpha = 1.0e-3_wp
ion_diss = Zero
flag_diss = .false.
cgamma = 300.0_wp/auToCm
HRSO = .false.
kext = .false.
! full density matrix saving time step
time_fdm = One/auToFs
! general idea is that additional features are disabled by default
! except for pulse flag
flag_so = .false.
flag_decay = .false.
flag_fdm = .false.
flag_dipole = .false.
flag_emiss = .false.
! number of incoming pulses N_pulse supposed to be 1 by default
! later when reading input it can be changed with
! reallocation of all corresponding arrays
flag_pulse = .true.
flag_acorrection = .false.
pulse_type = 'GAUSS'
N_pulse = 1
call mma_allocate(amp,N_pulse)
call mma_allocate(taushift,N_pulse)
call mma_allocate(pulse_vector,N_pulse,3)
call mma_allocate(sigma,N_pulse)
call mma_allocate(omega,N_pulse)
call mma_allocate(phi,N_pulse)
amp(1) = 2.5_wp
taushift(1) = Three/auToFs
pulse_vector(1,1) = cOne
pulse_vector(1,2) = cZero
pulse_vector(1,3) = cZero
sigma(1) = One/auToFs
omega(1) = Ten/autoev
phi(1) = Zero
linear_chirp = Zero

end subroutine rhodyn_init
