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
! Copyright (C) 2021-2023, Vladislav Kochetov                          *
!***********************************************************************

module rhodyn_data
!***********************************************************************
! Purpose: declaration of variables used between different subroutines *
!***********************************************************************
! i_rasscf         integer key used inside program (private):
!                  = 0 - nothing has been read from RASSCF
!                  = 1 - from RASSCF the spin-free H obtained
!                  = 2 - from RASSCF the CIs and Es obtained
!                  = 3 - from CASPT2 CIs and Es
!  runmode         integer key in input file:
!                  = 1 - start from rasscf/rassi input to propagation
!                  = 2 - start with sdprep.h5 file to propagation
!                  = 3 - just prepare sdprep.h5 file without propagation
!                  = 4 - charge migration (only rassi input)
!  N               nr of spin manifolds
!  Nstate          nr of total SF states acccounting for spin-degeneracy
!  ndet            dim(N) nr of the DET   for different spin manifolds
!  nconf           dim(N) nr of the CSFs for different spin manifolds
!  lroots          dim(N) nr of states (roots) for different spin manifolds
!  ispin           dim(N) multiplicity for different spin manifolds
!  istates         set of states (dim(Nstate)) included in dynamics
!  NDET_TOT        total nr of DETs (NDET_TOT=sum(NDET(I)))
!  nconftot        total nr of CSFs (nconftot=sum(NCONF(I)))
!  lrootstot       total nr of states (LROOTSTOT=sum(LROOTS(I)*ISPIN(I))
!  d               dimension of propagated matrices
!  N_populated     nr of the CSF/DET/SF/SO states to be populated
!  CI              CI(nconf,lroots,N) CI coefficients matrix including
!                      all spin manifolds
!  U_CI            U_CI(nconftot,lrootstot) the transformation matrix
!                      from spin-free states accounting for
!                       the spin-degeneracy  to CSFs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CSF ->  SF     !  U_CI^T     * . * U_CI                             !
!  SF  ->  CSF    !  U_CI       * . * U_CI^T                           !
!  SF  ->  SO     !  SO_CI^C    * . * SO_CI                            !
!  SO  ->  SF     !  SO_CI      * . * SO_CI^C                          !
!  SO  ->  CSF    !  CSF2SO     * . * CSF2SO^C                         !
!  CSF ->  SO     !  CSF2SO^C   * . * CSF2SO                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DTOC            DTOC(NDET,NCONF,N) transformation matrix
!                      from CSF to DET for particular spin manifold
!  DET2CSF         DET2CSF(NDET_TOT,NCONFTOT) full transform matrix
!                      from CSF to DET
!  V_SO            V_SO(NSS,NSS) SO-Hamiltonian matrix elements over
!                      spin components of spin-free eigenstates (SFS)
!  V_CSF           V_CSF(NCONFTOT,NCONFTOT) SO-Hamiltonian over CSFs
!  HTOTRE_CSF      real matrix, Hamiltonian in CSF basis of all spin
!                      manifolds which is extracted from RASSCF output
!  HTOT_CSF        complex matrix, HTOT_CSF=HTOTRE_CSF+V_CSF, total
!                       Hamiltonian in CSF basis including SO coupling
!  HSOCX           complex SO-Hamiltonian (taken from rassi)
!  p_style         what's populated: CSF, DET, SF or SO eigenstates
!  dysamp          matrix of Dyson amplitudes read from rassi (SO)
!  dysamp_bas      same as above in required basis
!
! out2_fmt, out3_fmt: formats for output writing declared in propagate
!
! k_max            max rank of spherical tensors used in propagation
! len_sph          length of spherical basis: different (k,q) pairs
!***********************************************************************

use Constants, only: kBoltzmann, auTokJ
use Definitions, only: wp, iwp

implicit none
private

! some abstract interfaces, which fit to several subroutines
! abstract interface
!   subroutine pulse_func(h0,ht,time,pcount)
!     import :: wp, iwp
!     complex(kind=wp), intent(in) :: h0(:,:)
!     complex(kind=wp), intent(out) :: ht(:,:)
!     real(kind=wp), intent(in) :: time
!     integer(kind=iwp), intent(in) :: pcount
!   end subroutine pulse_func
! end interface

! list of constants
real(kind=wp), parameter :: k_B = kBoltzmann/(auTokJ*1.0e3_wp), threshold = 1.0e-6_wp
integer(kind=iwp) :: ipglob, i_rasscf, runmode, N, Nstate, d, n_freq, ndet_tot, nconftot, lrootstot, maxnconf, maxlroots, n_sf
integer(kind=iwp), allocatable :: ispin(:), istates(:), lroots(:), nconf(:), ndet(:)
real(kind=wp), allocatable :: a_einstein(:,:), CI(:,:,:), DTOC(:,:,:), dysamp(:,:), E(:,:), H_CSF(:,:,:), HTOTRE_CSF(:,:), U_CI(:,:)
complex(kind=wp), allocatable :: CSF2SO(:,:), dipole(:,:,:), dipole_basis(:,:,:), DM0(:,:), dysamp_bas(:,:), E_SF(:), E_SO(:), &
                                 HSOCX(:,:), HTOT_CSF(:,:), SO_CI(:,:), tmp(:,:), U_CI_compl(:,:), V_CSF(:,:), V_SO(:,:)
! ----------------------------------------------------------------------
logical(kind=iwp) :: flag_decay, flag_dipole, flag_diss, flag_dyson, flag_emiss, flag_fdm, flag_so
character(len=256) :: basis, dm_basis, method, p_style, pulse_type
character(len=64) :: out2_fmt, out3_fmt, out_fmt, out_fmt_csf
character(len=6), allocatable :: rassd_list(:)
character(len=*), parameter :: int2real = '(a,2i5,2f16.8)', &
                               scha = '(1x,a,t52,a)', &
                               scmp = '(1x,a,t45,f5.2,sp,f5.2,"i")', &
                               sdbl = '(1x,a,t45,f9.3)', &
                               sint = '(1x,a,t45,i8)', &
                               slog = '(1x,a,t45,l8)'
! ----------------------------------------------------------------------
! prep, out hdf5 files (description in cre_prep, cre_out)
integer(kind=iwp) :: out_decay_i, out_decay_r, out_dm_csf, out_dm_sf, out_dm_so, out_emiss, out_fdmi, out_fdmr, out_freq, &
                     out_ham_i, out_ham_r, out_id, out_pulse, out_t, out_tfdm, out_tout, prep_ci, prep_csfsoi, prep_csfsor, &
                     prep_dipolei, prep_dipoler, prep_dm_i, prep_dm_r, prep_do, prep_fhi, prep_fhr, prep_hcsf, prep_id, prep_uci, &
                     prep_vcsfi, prep_vcsfr
! predefined units for writing output, think of more clever choice:
integer(kind=iwp) :: lu_csf = 32, &
                     lu_dip = 34, &
                     lu_pls = 33, &
                     lu_sf = 31, &
                     lu_so = 30
! ----------------------------------------------------------------------
! variables used for density matrix propagation
integer(kind=iwp) :: N_L2, N_L3, N_Populated, Nmode, Npop, Nstep, Ntime_tmp_dm, Nval
logical(kind=iwp) :: HRSO, kext
real(kind=wp) :: alpha, dt, errorthreshold, finaltime, cgamma, initialtime, ion_diss, safety, T, tau_L2, tau_L3, time_fdm, &
                 timestep, tout
real(kind=wp), allocatable :: dgl(:), emiss(:)
complex(kind=wp), allocatable :: decay(:,:), density0(:,:), densityt(:,:), hamiltonian(:,:), hamiltoniant(:,:), k_bar_basis(:,:), &
                                 kab_basis(:,:), pulse_vector(:,:)
! Runge-Kutta midpoints
complex(kind=wp), allocatable :: ak1(:,:), ak2(:,:), ak3(:,:), ak4(:,:), ak5(:,:), ak6(:,:)
! pulse characteristics
integer(kind=iwp) :: N_pulse, power_shape
logical(kind=iwp) :: flag_acorrection, flag_pulse
real(kind=wp) :: linear_chirp
complex(kind=wp) :: pulse_vec(3)
real(kind=wp), allocatable :: amp(:), omega(:), phi(:), sigma(:), taushift(:)
! ----------------------------------------------------------------------
! for the propagation of rho in spherical tensors basis
integer(kind=iwp) :: k_max, len_sph, q_max
integer(kind=iwp), allocatable :: k_ranks(:), list_sf_mult(:), list_sf_states(:), list_so_mult(:), list_so_sf(:), q_proj(:)
real(kind=wp), allocatable :: list_sf_spin(:), list_so_proj(:), list_so_spin(:)
complex(kind=wp), allocatable :: mirr(:,:)
complex(kind=wp), allocatable :: midk1(:,:,:), midk2(:,:,:), midk3(:,:,:), midk4(:,:,:), V_SO_red(:,:,:)
complex(kind=wp), allocatable :: Y1(:,:,:), Y2(:,:,:)

public :: a_einstein, ak1, ak2, ak3, ak4, ak5, ak6, alpha, amp, basis, cgamma, CI, CSF2SO, d, decay, density0, densityt, dgl, &
          dipole, dipole_basis, DM0, DM_basis, dt, DTOC, dysamp, dysamp_bas, E, E_SF, E_SO, emiss, errorthreshold, finaltime, &
          flag_acorrection, flag_decay, flag_dipole, flag_diss, flag_dyson, flag_emiss, flag_fdm, flag_pulse, flag_so, H_CSF, &
          hamiltonian, hamiltoniant, HRSO, HSOCX, HTOT_CSF, HTOTRE_CSF, i_rasscf, initialtime, int2real, ion_diss, ipglob, ispin, &
          istates, k_b, K_bar_basis, k_max, k_ranks, kab_basis, kext, len_sph, linear_chirp, list_sf_mult, list_sf_spin, &
          list_sf_states, list_so_mult, list_so_proj, list_so_sf, list_so_spin, lroots, lrootstot, lu_csf, lu_dip, lu_pls, lu_sf, &
          lu_so, maxlroots, maxnconf, method, midk1, midk2, midk3, midk4, mirr, N, n_freq, N_L2, N_L3, N_Populated, N_pulse, n_sf, &
          nconf, nconftot, ndet, ndet_tot, Nmode, Npop, Nstate, Nstep, Ntime_tmp_dm, Nval, omega, out2_fmt, out3_fmt, out_decay_i, &
          out_decay_r, out_dm_csf, out_dm_sf, out_dm_so, out_emiss, out_fdmi, out_fdmr, out_fmt, out_fmt_csf, out_freq, out_ham_i, &
          out_ham_r, out_id, out_pulse, out_t, out_tfdm, out_tout, p_style, phi, power_shape, prep_ci, prep_csfsoi, prep_csfsor, &
          prep_dipolei, prep_dipoler, prep_dm_i, prep_dm_r, prep_do, prep_fhi, prep_fhr, prep_hcsf, prep_id, prep_uci, prep_vcsfi, &
          prep_vcsfr, pulse_type, pulse_vec, pulse_vector, q_max, q_proj, rassd_list, runmode, safety, scha, scmp, sdbl, sigma, &
          sint, slog, SO_CI, T, tau_L2, tau_L3, taushift, threshold, time_fdm, timestep, tmp, tout, U_CI, U_CI_compl, V_CSF, V_SO, &
          V_SO_red, Y1, Y2

end module rhodyn_data
