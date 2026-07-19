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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!               2026, Bhavnesh Jangid                                  *
!***********************************************************************

module lpdft
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: wp, iwp, u6
  use general_data, only: ispin, jobiph, jobold, nash, nbas, nfro, nish, nsym, ntot1, ntot2
  use mcpdft_input, only: mcpdft_options
  use rasscf_global, only: DFTFOCK, ExFac, IADR15, NACPAR, NACPR2, NFINT, NonEq, PotNuc, NROOTS, WEIGHT
  use rctfld_module, only: lRF
  implicit none
  private
  logical :: do_lpdft=.false.

  real(kind=wp) :: hconst = 0.0D0
  real(kind=wp) :: e_ot_lpdft = 0.0D0
  real(kind=wp), allocatable, dimension(:) :: dm1_core
  real(kind=wp), allocatable, dimension(:) :: casd1_0, casd1s_0, casd2_0
  real(kind=wp), allocatable, dimension(:) :: veff1, veff2
  real(kind=wp), allocatable, dimension(:) :: h1e_lpdft_ao, h1e_lpdft_mo, h2e_lpdft_mo

  public :: do_lpdft
  public :: lpdft_kernel
  contains

  ! LPDFT Kernel Function.
  ! Basically, compute the L-PDFT Hamiltonian, and then contract it with transition density
  ! matrices to compute the effective Hamiltonian, whose eigen vectors are state-interaction vectors
  ! and will be used to roate the SA-RAS CI vectors.

  subroutine lpdft_kernel(CMO)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    call init_vars(CMO)

    call set_iadr15()
    ! Want to first get the zero-order densities
    call get_zero_order_densities(casd1_0, casd1s_0, casd2_0)
    call compute_Eot_and_potentials(CMO)
    call compute_hconst(CMO)
    call compute_h1e_lpdft()
    call compute_h2e_lpdft()

    ! Using hconst, h1e, and h2e now we can compute the effective Hamiltonian and then diagonalize
    ! it to get the state-interaction vectors
    call release_vars()
  end subroutine lpdft_kernel

    subroutine get_zero_order_densities(casd1_0, casd1s_0, casd2_0)
      real(kind=wp), dimension(nacpar), intent(out) :: casd1_0, casd1s_0
      real(kind=wp), dimension(nacpr2), intent(out) :: casd2_0

      real(kind=wp), dimension(nacpar) :: casd1, casd1s
      real(kind=wp), dimension(nacpr2) :: casd2
      integer(kind=iwp) :: disk, root

      disk = iAdr15(3)
      casd1_0(:) = 0.0D0
      casd1s_0(:) = 0.0D0
      casd2_0(:) = 0.0D0

      do root=1, nroots
        call dDaFile(JOBOLD,2, casd1, nacpar, disk)
        call dDaFile(JOBOLD,2, casd1s, nacpar, disk)
        call dDaFile(JOBOLD,2, casd2, nacpr2, disk)
        ! Need a dummy read for some reason
        call dDaFile(JOBOLD,0, casd2, nacpr2, disk)

        call daxpy_(nacpar, weight(root), casd1, 1, casd1_0, 1)
        call daxpy_(nacpar, weight(root), casd1s, 1, casd1s_0, 1)
        call daxpy_(nacpr2, weight(root), casd2, 1, casd2_0, 1)
      end do
    end subroutine get_zero_order_densities

    subroutine compute_Eot_and_potentials(CMO)
      use KSDFT_Info, only: do_pdftpot

      real(kind=wp), dimension(ntot2), intent(in) :: CMO

      ! Some dummy vars that don't matter
      real(kind=wp), dimension(ntot1) :: htmp, gtmp
      ! The density matrices
      real(kind=wp), dimension(ntot2) :: casd1AO_0, casd1sAO_0
      real(kind=wp), dimension(ntot1) :: d1ao_0, casD1AO_folded_0, d1sao_0
      integer(kind=iwp) :: charge
      integer(kind=iwp), external :: get_charge

      logical, parameter :: First=.True., Dff=.False., Do_DFT=.True.

      ! Convert CASD1 (MO) -> AO Basis
      call Get_D1A_RASSCF(CMO, casd1_0, casd1AO_0)
      call Get_D1A_RASSCF(CMO, casd1s_0, casd1sAO_0)

      call Fold(nSym, nBas, dm1_core, d1ao_0)
      call Fold(nSym, nBas, casD1AO_0, casD1AO_folded_0)
      call Fold(nSym, nBas, casd1sAO_0, d1sao_0)
      ! Create full 1RDM (inactive/frozen + active) and store in d1ao_0
      call daxpy_(nTot1, 1.0D0, casD1AO_folded_0, 1, d1ao_0, 1)

      ! Load in the nuclear potential
      Call Get_dScalar('PotNuc',potNuc)

      charge = get_charge()

      ! We need the 1 and 2-electron potentials
      do_pdftpot = .True.

      veff1 = 0.0D0
      veff2 = 0.0D0
      call put_dArray('ONTOPO', veff1, ntot1)
      call put_dArray('ONTOPT', veff2, nfint)
      call put_dArray('D1ao', d1ao_0, ntot1)
      call put_dArray('D1Sao', d1sao_0, ntot1)
      call put_dArray('D1mo', casd1_0, nacpar)
      call put_dArray('P2mo', casd2_0, nacpr2)
      call put_iArray('nFro', nfro, nsym)
      call put_iArray('nAsh', nash, nsym)
      call put_iArray('nIsh', nish, nsym)

      htmp = 0.0D0
      gtmp = 0.0D0
      call DrvXV(htmp, gtmp, d1ao_0, potnuc, nTot1, First, Dff, NonEq, &
                 lRF, mcpdft_options%otfnal%otxc, ExFac, charge, ispin, DFTFOCK, Do_DFT)

      ! Load in the on-top potentials
      call get_dArray("ONTOPO", veff1, ntot1)
      call get_dArray("ONTOPT", veff2, nfint)

      ! Read the on-top energy
      call get_dScalar('CASDFT energy', e_ot_lpdft)
    end subroutine compute_Eot_and_potentials

  subroutine compute_hconst(CMO)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    real(kind=wp) :: e1_lpdft, e2_lpdft, eCoulLPDFT, eVeffLPDFT, hnuc
    real(kind=wp), allocatable, dimension(:) :: casd1AO_0, casD1AO_folded_0, d1ao_0, veff2_act
    real(kind=wp), allocatable, dimension(:) :: veff_poten, veff_poten_act, vj
    real(kind=wp), external :: dDot_

    ! Focusing only on the non-hybrid for the time-being.

    ! hconst = hnuc + E_OT -0.5 * g_pqrs D_pq D_rs - V_pq D_pq -1/2 * V_pqrs D_pqrs

    hnuc = 0.0D0
    e1_lpdft = 0.0D0
    e2_lpdft = 0.0D0
    eCoulLPDFT = 0.0D0
    eVeffLPDFT = 0.0D0

    ! Get the nuclear repulsion energy
    Call Get_dScalar('PotNuc',hnuc)

    call mma_allocate(casd1AO_0, ntot2, "CASD1AO_0")
    call mma_allocate(casD1AO_folded_0, ntot1, "CASD1AO_FOLDED_0")
    call mma_allocate(d1ao_0, ntot1, "D1AO_0")
    call mma_allocate(vj, ntot1, "VJ_LPDFT")
    call mma_allocate(veff_poten, ntot1, "VEFF_POTEN_LPDFT")
    call mma_allocate(veff_poten_act, nacpar, "VEFF_POTEN_ACT_LPDFT")
    call mma_allocate(veff2_act, nacpr2, "VEFF2_ACT")

    ! Compute the 1e contribution
    call Get_D1A_RASSCF(CMO, casd1_0, casd1AO_0)
    call Fold(nSym, nBas, dm1_core, d1ao_0)
    call Fold(nSym, nBas, casd1AO_0, casD1AO_folded_0)
    call daxpy_(ntot1, 1.0D0, casD1AO_folded_0, 1, d1ao_0, 1)

    call get_coulomb(CMO, dm1_core, casd1AO_0, vj)
    call get_veff_poten(CMO, dm1_core, casd1AO_0, veff_poten)

    eCoulLPDFT = 0.5D0*dDot_(ntot1, vj, 1, d1ao_0, 1)
    e1_lpdft = dDot_(ntot1, veff1, 1, d1ao_0, 1)

    call ao2moact_1e(CMO, veff_poten, veff_poten_act, nSym, nBas, nAsh, nFro, nIsh)

    eVeffLPDFT = 2.0D0*dDot_(nacpar, veff_poten_act, 1, casd1_0, 1)
    e1_lpdft = e1_lpdft+eVeffLPDFT

    veff2_act(:) = 0.0D0
    call Get_TUVX(veff2, veff2_act)

    e2_lpdft = 0.5D0*dDot_(nacpr2, veff2_act, 1, casd2_0, 1)

    hconst = hnuc + e_ot_lpdft - eCoulLPDFT - e1_lpdft - e2_lpdft

    write(u6,'(6X,A,ES20.10)') 'LPDFT nuclear repulsion energy: ', hnuc
    write(u6,'(6X,A,ES20.10)') 'LPDFT on-top energy:           ', e_ot_lpdft
    write(u6,'(6X,A,ES20.10)') 'LPDFT Coulomb double count:    ', eCoulLPDFT
    write(u6,'(6X,A,ES20.10)') 'LPDFT effective potential term:', eVeffLPDFT
    write(u6,'(6X,A,ES20.10)') 'LPDFT one-body OT double count:', e1_lpdft
    write(u6,'(6X,A,ES20.10)') 'LPDFT two-body OT double count:', e2_lpdft
    write(u6,'(6X,A,ES20.10)') 'LPDFT scalar Hamiltonian term: ', hconst

    call mma_deallocate(veff2_act)
    call mma_deallocate(veff_poten_act)
    call mma_deallocate(veff_poten)
    call mma_deallocate(vj)
    call mma_deallocate(d1ao_0)
    call mma_deallocate(casD1AO_folded_0)
    call mma_deallocate(casd1AO_0)
  end subroutine compute_hconst

  subroutine compute_h1e_lpdft()
    ! Compute the hcore:
    ! Add the vot potentials to the hcore:
    ! Transform the hcore to the MO basis:
  end subroutine compute_h1e_lpdft

  subroutine compute_h2e_lpdft()
    ! Compute the h2e:
    ! Transform the h2e to the MO basis:
  end subroutine compute_h2e_lpdft

  subroutine set_iadr15()
    ! Set the iAdr15 variable in rasscf.fh so I don't have to do any
    ! weird jobiph/jobold looking.
    logical :: found
    integer :: iad15

    iadr15 = 0
    iad15 = 0
    call f_inquire('JOBOLD', Found)
    if(.not. found) then
      call f_inquire('JOBIPH', Found)
      if(found) then
        JOBOLD=JOBIPH
      end if
    end if
    if(found) then
      if(JOBOLD <= 0) then
        JOBOLD = 20
        call DaName(JOBOLD, 'JOBOLD')
      end if
    end if
    call iDaFile(JOBOLD, 2, iadr15, 15, iad15)
  end subroutine set_iadr15

  subroutine init_vars(CMO)
    real(kind=wp), dimension(ntot2) :: CMO

    call mma_allocate(casd1_0, nacpar, "CASD1_0")
    call mma_allocate(casd1s_0, nacpar, "CASD1s_0")
    call mma_allocate(casd2_0, nacpr2, "CASD2_0")
    call mma_allocate(dm1_core, ntot2, "dm1_core")
    call mma_allocate(veff1, ntot1, "1E-OnTopPot")
    call mma_allocate(veff2, nfint, "2E-OnTopPot")
    call mma_allocate(h1e_lpdft_ao, ntot1, "H1E_LPDFT_AO")
    call mma_allocate(h1e_lpdft_mo, nacpar, "H1E_LPDFT_MO")
    call mma_allocate(h2e_lpdft_mo, nacpr2, "H2E_LPDFT_MO")

    write(u6,*) "size(veff1) =", size(veff1)
    write(u6,*) "size(veff2) =", size(veff2)
    write(u6,*) "size(casdm1_0) =", size(casd1_0)
    write(u6,*) "size(casdm1s_0) =", size(casd1s_0)
    write(u6,*) "size(casdm2) =", size(casd2_0)
    write(u6,*) "size(dm1_core) =", size(dm1_core)
    write(u6,*) "size(h1e_lpdft_ao) =", size(h1e_lpdft_ao)
    write(u6,*) "size(h1e_lpdft_mo) =", size(h1e_lpdft_mo)
    write(u6,*) "size(h2e_lpdft_mo) =", size(h2e_lpdft_mo)

    call Get_D1I_RASSCF(CMO, dm1_core)

  end subroutine init_vars

  subroutine release_vars()
    call mma_deallocate(casd1_0)
    call mma_deallocate(casd1s_0)
    call mma_deallocate(casd2_0)
    call mma_deallocate(h1e_lpdft_ao)
    call mma_deallocate(h1e_lpdft_mo)
    call mma_deallocate(h2e_lpdft_mo)
    call mma_deallocate(dm1_core)
    call mma_deallocate(veff1)
    call mma_deallocate(veff2)
  end subroutine release_vars

end module lpdft
