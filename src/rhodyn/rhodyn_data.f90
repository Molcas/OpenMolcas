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
module rhodyn_data
  use definitions, only: wp, iwp
  use constants, only: kBoltzmann, auTokJ
  implicit none
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

! out_fmt,out1_fmt:formats for output writing declared in propagate.f
!***********************************************************************
! some abstract interfaces, which fit to several subroutines
  abstract interface
    subroutine rk_fixed_step(t0,y)
      real(8), intent(in) :: t0
      complex(8), dimension(:,:), intent(inout) :: y
    end subroutine
    subroutine rk_adapt_step(t0,y,err)
      real(8), intent(in) :: t0
      real(8), intent(out) :: err
      complex(8), dimension(:,:), intent(inout) :: y
    end subroutine
    subroutine pulse_func(h0,ht,time,count)
      complex(8), dimension(:,:), intent(in) :: h0
      complex(8), dimension(:,:), intent(out) ::ht
      real(8), intent(in) :: time
      integer, intent(in), optional :: count
    end subroutine pulse_func
    subroutine equation_func(time,rho_t,res)
      real(8), intent(in) :: time
      complex(8), dimension(:,:), intent(in) :: rho_t
      complex(8), dimension(:,:), intent(out) :: res
    end subroutine
  end interface
  ! list of dummy integers
  integer(kind=iwp) :: i,j,k,l,ii,jj,kk,ll
  ! list of constants
  real(kind=wp), parameter :: k_B = kBoltzmann/(auTokJ*1.0e3_wp) ,&
                              threshold = 1.0d-06                ,&
                              tiny      = 1.0d-20
!                        k_B  = 3.1668114d-6 Hartree/K
!                        Debyetoau = 0.393456
!                        autoev = 27.211396132
!                        cmtoau    = 4.5563d-6
!                        fstoau    = 41.3393964d0
  complex(kind=wp),parameter:: zero = (0.0d0,0.0d0)    ,&
                               one  = (1.0d0,0.0d0)    ,&
                               onei = (0.0d0,1.0d0)
  integer(kind=iwp) :: ireturn, ipglob, error, i_rasscf ,&
                       runmode ,N, Nstate, d, n_freq,&
                       ndet_tot, nconftot, lrootstot    ,&
                       maxnum, maxnconf, maxlroots
  integer(kind=iwp), dimension(:), allocatable :: ndet,nconf,lroots,&
                                                  ispin,istates
  real(kind=wp), dimension(:,:,:), allocatable :: H_CSF,CI,DTOC
  real(kind=wp), dimension(:,:),  allocatable :: E, U_CI, HTOTRE_CSF,&
                                                 dysamp,a_einstein
  complex(kind=wp), dimension(:,:,:),allocatable:: dipole,dipole_basis
  complex(kind=wp), dimension(:,:),  allocatable:: V_SO,V_CSF,tmp,&
                                              HTOT_CSF,CSF2SO,DM0,&
                                              U_SO,SO_CI,HSOCX,&
                                              U_CI_compl,dysamp_bas
  complex(kind=wp), dimension(:), allocatable :: E_SF, E_SO
! ---------------------------------------------------------------------
  logical :: flag_so, flag_decay, flag_diss, flag_fdm, &
             flag_dyson, flag_emiss, flag_test, flag_dipole
  character(len=256) :: pulse_type, dm_basis, p_style, method, basis
  character(len=6),dimension(:),allocatable:: rassd_list, hr_list
  character(len=*),parameter :: sint='(x,a,t45,i8)'              ,&
                                scha='(x,a,t52,a)'               ,&
                                sdbl='(x,a,t45,f9.3)'            ,&
                                scmp='(x,a,t45,f5.2,sp,f5.2,"i")',&
                                slog='(x,a,t45,l8)'              ,&
                                int2real='(a,2i5,2f16.8)'
  character(len=64) :: out_fmt, out1_fmt, out_fmt_csf, out1_fmt_csf,&
                       out2_fmt, out3_fmt
! ---------------------------------------------------------------------
! prep, out hdf5 files (description in cre_prep.f, cre_out.f)
  integer(kind=iwp) :: prep_id, &
                    prep_ci,  prep_dipoler, prep_dm_r,  &
                    prep_uci, prep_dipolei, prep_dm_i,  &
                    prep_utu, prep_fullh,   prep_vcsfr, &
                    prep_fhr, prep_csfsor,  prep_vcsfi, &
                    prep_fhi, prep_csfsoi,  prep_hcsf,  &
                    prep_do,                                &
                    out_id,  &
                    out_dm_so, out_dm_sf, out_dm_csf,   &
                    out_pulse, out_tout, out_t, out_fdm,&
                    out_decay_r, out_decay_i, out_ham_r,&
                    out_ham_i, out_freq, out_emiss,     &
                    out_tfdm
! predefined units for writing output, think of more clever choice:
  integer(kind=iwp):: lu_so =30, &
                      lu_sf =31, &
                      lu_csf=32, &
                      lu_pls=33, &
                      lu_dip=34
! ---------------------------------------------------------------------
! variables used for density matrix propagation
  integer(kind=iwp)    :: N_Populated,Ntime_tmp_dm,Nstep,Npop
  integer(kind=iwp)    :: Nval,N_L3,N_L2,Nmode
  logical              :: HRSO, kext
  logical,dimension(5) :: ion_blocks ! for processing Dyson matrix
  real(kind=wp)        :: T,tau_L3,tau_L2,gamma,&
              initialtime,finaltime,timestep,dt,deltaE,V,tout,&
              time_fdm, errorthreshold, alpha, safety, ion_diss
  real(kind=wp),dimension(:),allocatable :: dgl, emiss
  real(kind=wp),dimension(6) :: temp_vec
  complex(8),dimension(:,:),allocatable :: decay, pulse_vector,&
                                 density0,densityt            ,&
                                 hamiltonian,hamiltoniant     ,&
                                 kab_basis, k_bar_basis
! Runge-Kutta midpoints
  complex(8),dimension(:,:),allocatable::ak1,ak2,ak3,ak4,ak5,ak6,y5
! pulse characteristics
  logical :: flag_pulse, flag_acorrection
  integer(kind=iwp) :: N_pulse, power_shape
  real(kind=wp) :: linear_chirp, t_local, omega_local
  real(kind=wp), dimension(:), allocatable :: amp, omega, sigma ,&
                                              phi, taushift
  complex(kind=wp), dimension(3) :: pulse_vec, E_field
end module rhodyn_data
