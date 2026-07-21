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
  use general_data, only: ispin, jobiph, jobold, nash, nbas, nfro, nish, norb, nsym, ntot1, ntot2
  use mcpdft_input, only: mcpdft_options
  use rasscf_global, only: dftfock, exfac, iadr15, istord, istorp, lroots, nac, nacpar, nacpr2, nfint, noneq, potnuc, nroots, &
                           weight
  use rctfld_module, only: lRF

  implicit none
  private
  logical :: do_lpdft=.false.

  real(kind=wp) :: hconst = 0.0D0
  real(kind=wp) :: e_ot_lpdft = 0.0D0
  real(kind=wp) :: e_nuc_lpdft = 0.0D0
  real(kind=wp) :: e_wfn_lpdft = 0.0D0
  real(kind=wp), allocatable, dimension(:) :: dm1_core
  real(kind=wp), allocatable, dimension(:) :: casd1_0, casd1s_0, casd2_0
  real(kind=wp), allocatable, dimension(:) :: veff1, veff2
  real(kind=wp), allocatable, dimension(:) :: h1e_lpdft_ao, h1e_lpdft_act, h2e_lpdft_act
  real(kind=wp), allocatable, dimension(:) :: vj, vjk

  public :: do_lpdft
  public :: lpdft_kernel
  contains

  !***********************************************************************
  !> LPDFT Kernel Function.
  !> Basically, compute the L-PDFT Hamiltonian, and then contract it with transition density
  !> matrices to compute the effective Hamiltonian, whose eigen vectors are state-interaction vectors
  !> and will be used to roate the SA-RAS CI vectors.
  !>
  !***********************************************************************

  subroutine lpdft_kernel(CMO,e_mcscf)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    real(kind=wp), intent(in) :: e_mcscf(:)

    call init_vars(CMO)

    call set_iadr15()

    call get_zero_order_densities(casd1_0, casd1s_0, casd2_0)
    call compute_Eot_and_potentials(CMO)

    call compute_hconst(CMO)
    call compute_h1e_lpdft(CMO)
    call compute_h2e_lpdft()

    ! Using hconst, h1e, and h2e now we can compute the effective Hamiltonian and then diagonalize
    ! it to get the state-interaction vectors
    call construct_lpdft_ham()

    call print_lpdft_finalize(e_nuc_lpdft,e_wfn_lpdft,e_ot_lpdft,e_mcscf,size(e_mcscf))

    call release_vars()
  end subroutine lpdft_kernel

  !***********************************************************************
  !> get_zero_order_densities
  !> Basically, the D0 = sum_i w_i D_i, where D_i is the density matrix
  !> for the i-th root and w_i is the weight of the i-th root.
  !>***********************************************************************
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

  !>***********************************************************************
  !> compute_Eot_and_potentials
  !> Compute the on-top energy and the on-top potentials. This is done by
  !> calling the DrvXV routine, which computes the effective potential for
  !> a given density matrix.
  !>***********************************************************************
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

      veff1(:) = 0.0D0
      veff2(:) = 0.0D0
      vjk(:) = 0.0D0
      call put_dArray('ONTOPO', veff1, ntot1)
      call put_dArray('ONTOPT', veff2, nfint)
      call put_dArray('FI_V', vjk, ntot1)
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
                 lRF, mcpdft_options%otfnal%otxc, ExFac, charge, ispin, dftfock, Do_DFT)

      ! Load in the on-top potentials.
      call get_dArray("ONTOPO", veff1, ntot1)
      call get_dArray("ONTOPT", veff2, nfint)
      call get_dArray("FI_V", vjk, ntot1)

      ! Read the on-top energy
      call get_dScalar('CASDFT energy', e_ot_lpdft)
  end subroutine compute_Eot_and_potentials


  !***********************************************************************
  !> Compute core density matrix in MO basis. Note: this would be a diagonal matrix.
  !> I need to write this wrapper function to make sure I am following the OpenMolcas
  !> convention of storing the density matrix in triangular form, along with
  !> the symmetry blocks.
  !***********************************************************************
  subroutine get_core_dm1_mo(dm1_core_mo)
    use Index_Functions, only: iTri, nTri_Elem
    real(kind=wp), intent(out) :: dm1_core_mo(ntot1)
    integer(kind=iwp) :: iCore, iMOff, iOrb, iSym

    dm1_core_mo(:) = 0.0D0
    iMOff = 1
    do iSym=1,nSym
      iOrb = nOrb(iSym)
      do iCore=1,nIsh(iSym)
        dm1_core_mo(iMOff+iTri(iCore,iCore)-1) = 2.0D0
      end do
      iMOff = iMOff+nTri_Elem(iOrb)
    end do
  end subroutine get_core_dm1_mo

  !***********************************************************************
  !> get_dm1_mo
  !> Get the 1e density matrix in MO basis.
  !> dm1_mo = dm1_core_mo + dm1_active_mo
  !***********************************************************************
  subroutine get_dm1_mo(d1mo)
    use Index_Functions, only: iTri, nTri_Elem

    real(kind=wp), intent(out) :: d1mo(ntot1)
    integer(kind=iwp) :: iAct, iActOff, iCore, iMOff, iMO, iOrb, iSym, jAct, jMO

    d1mo(:) = 0.0D0
    iMOff = 1
    iActOff = 1
    do iSym=1,nSym
      iOrb = nOrb(iSym)
      do iCore=1,nIsh(iSym)
        d1mo(iMOff+iTri(iCore,iCore)-1) = 2.0D0
      end do
      do iAct=1,nAsh(iSym)
        iMO = nIsh(iSym)+iAct
        do jAct=1,iAct
          jMO = nIsh(iSym)+jAct
          d1mo(iMOff+iTri(iMO,jMO)-1) = casd1_0(iActOff)
          iActOff = iActOff+1
        end do
      end do
      iMOff = iMOff+nTri_Elem(iOrb)
    end do
  end subroutine get_dm1_mo

  !***********************************************************************
  !> get_e_vhf_core
  !> Get the core contribution to the vhf_c term.
  !> e_vhf_c_core = sum_{iCore} vhf_c(iCore,iCore)
  !> vhf_c is in MO basis, and is stored in triangular form.
  !***********************************************************************
  function get_e_vhf_core(vhf_c_mo) result(e_core)
    use Index_Functions, only: iTri, nTri_Elem

    real(kind=wp), intent(in) :: vhf_c_mo(ntot1)
    real(kind=wp) :: e_core
    integer(kind=iwp) :: iCore, iMOff, iOrb, iSym

    e_core = 0.0D0
    iMOff = 1
    do iSym=1,nSym
      iOrb = nOrb(iSym)
      do iCore=1,nIsh(iSym)
        e_core = e_core + vhf_c_mo(iMOff+iTri(iCore,iCore)-1)
      end do
      iMOff = iMOff+nTri_Elem(iOrb)
    end do
  end function get_e_vhf_core

  !***********************************************************************
  !> get_core_energy
  !> Get the core energy contribution to the L-PDFT Hamiltonian.
  !> Basically, the E_MCSCF = E_core + E_active.
  !> In this function, we compute the E_core term, which is given by:
  !> E_core = Tr(hcore*dm1_core) + Tr(veff1*dm1_core) + Tr(vj*dm1_core)
  !***********************************************************************
  function get_core_energy(CMO) result(e_core)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    real(kind=wp) :: e_core
    real(kind=wp), allocatable, dimension(:) :: dm1_core_mo, hcore_mo, vj_mo
    real(kind=wp), external :: dDot_

    call mma_allocate(dm1_core_mo, ntot1, "DM1_CORE_MO")
    call mma_allocate(hcore_mo, ntot1, "HCORE_MO_LPDFT")
    call mma_allocate(vj_mo, ntot1, "VJ_MO_LPDFT")

    e_core = 0.0D0
    dm1_core_mo(:) = 0.0D0
    hcore_mo(:) = 0.0D0
    vj_mo(:) = 0.0D0
    h1e_lpdft_ao(:) = 0.0D0

    call get_core_dm1_mo(dm1_core_mo)
    call get_hcore(h1e_lpdft_ao)
    call ao2mo_1e(CMO, h1e_lpdft_ao, hcore_mo, nSym, nBas, nOrb, nFro)
    call ao2mo_1e(CMO, vj, vj_mo, nSym, nBas, nOrb, nFro)

    e_core = dDot_(ntot1, hcore_mo, 1, dm1_core_mo, 1)
    e_core = e_core + dDot_(ntot1, veff1, 1, dm1_core_mo, 1)
    e_core = e_core + dDot_(ntot1, vj_mo, 1, dm1_core_mo, 1)

    call mma_deallocate(vj_mo)
    call mma_deallocate(hcore_mo)
    call mma_deallocate(dm1_core_mo)
  end function get_core_energy


  !>***********************************************************************
  !> compute_hconst:
  !> the energy contribution from the weighted average densities.
  !> This is done by computing the nuclear repulsion energy, the core energy,
  !> the on-top energy, and the one- and two-electron contributions from the
  !> on-top potentials. The constant term is then given by:
  !> hconst = E_nuc + E_core + E_ot - E_1 - E_2
  !>***********************************************************************
  subroutine compute_hconst(CMO)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    real(kind=wp) :: e1_lpdft, e2_lpdft, e_core_lpdft, e_coul_lpdft, e_vhf_c, e_vhf_c_core, hnuc
    real(kind=wp), allocatable, dimension(:) :: casd1AO_0, casD1AO_folded_0, d1ao_0, d1mo_0
    real(kind=wp), allocatable, dimension(:) :: vjk_act
    real(kind=wp), external :: dDot_

    ! Focusing only on the non-hybrid for the time-being.

    hnuc = 0.0D0
    e1_lpdft = 0.0D0
    e2_lpdft = 0.0D0
    e_core_lpdft = 0.0D0
    e_coul_lpdft = 0.0D0
    e_vhf_c = 0.0D0
    e_vhf_c_core = 0.0D0

    ! Term-1: Get the nuclear repulsion energy
    Call Get_dScalar('PotNuc',hnuc)

    call mma_allocate(casd1AO_0, ntot2, "CASD1AO_0")
    call mma_allocate(casD1AO_folded_0, ntot1, "CASD1AO_FOLDED_0")
    call mma_allocate(d1ao_0, ntot1, "D1AO_0")
    call mma_allocate(d1mo_0, ntot1, "D1MO_0")
    call mma_allocate(vjk_act, nacpar, "VJK_ACT_LPDFT")

    ! First compute the hconst energy then compute the core energy from the
    ! L-PDFT Hamiltonian, and finally adding these would be used as the constant
    ! energy term to shift the diagonal elements of the effective Hamiltonian.

    ! Core density matrix in AO basis
    call Get_D1A_RASSCF(CMO, casd1_0, casd1AO_0)
    call Fold(nSym, nBas, dm1_core, d1ao_0)
    ! Just the active space density matrix in AO basis
    call Fold(nSym, nBas, casd1AO_0, casD1AO_folded_0)
    ! Complete 1e density matrix in AO basis
    call daxpy_(ntot1, 1.0D0, casD1AO_folded_0, 1, d1ao_0, 1)

    call get_dm1_mo(d1mo_0)

    ! Compute the classical Coulomb energy
    vj(:) = 0.0D0
    call get_coulomb(CMO, dm1_core, casd1AO_0, vj)
    e_coul_lpdft = 0.5D0*dDot_(ntot1, vj, 1, d1ao_0, 1)

    ! Compute the veff1 (1e on-top potential) contribution
    e1_lpdft = dDot_(ntot1, veff1, 1, d1mo_0, 1)

    ! Finally add the 1e energy contr.
    e1_lpdft = e1_lpdft + e_coul_lpdft

    ! Compute the 2e contribution
    e2_lpdft = contract_puvx_2rdm(veff2, casd2_0)

    ! Compute the vhf_c contribution from active and core densities
    call mo2act_1e(vjk, vjk_act, nSym, nOrb, nAsh, nIsh)
    e_vhf_c = dDot_(nacpar, vjk_act, 1, casd1_0, 1)
    e_vhf_c_core = get_e_vhf_core(vjk)

    e2_lpdft = e2_lpdft + e_vhf_c + e_vhf_c_core

    ! Now computing the core-energy contributions. Keep the veff2 core term
    ! explicit here so it cancels the same term in e2_lpdft, matching the
    ! PySCF get_lpdft_hconst/transformed_h1e_for_cas split.
    e_core_lpdft = get_core_energy(CMO) + e_vhf_c_core

    ! Finally, compute the hconst term
    hconst = hnuc + e_ot_lpdft - e1_lpdft - e2_lpdft + e_core_lpdft
    e_nuc_lpdft = hnuc
    e_wfn_lpdft = hnuc + e_core_lpdft

    call mma_deallocate(vjk_act)
    call mma_deallocate(d1mo_0)
    call mma_deallocate(d1ao_0)
    call mma_deallocate(casD1AO_folded_0)
    call mma_deallocate(casd1AO_0)
  end subroutine compute_hconst

  function contract_puvx_2rdm(puvx, p2) result(energy)
    real(kind=wp), intent(in) :: puvx(nfint), p2(nacpr2)
    real(kind=wp) :: energy
    real(kind=wp), allocatable :: p2_reordered(:)
    integer(kind=iwp) :: isym, lrow, nao, nat, nio, no, nuvx, p_idx, puvx_idx

    energy = 0.0D0
    if (istorp(nsym+1) <= 0) return

    call mma_allocate(p2_reordered, istorp(nsym+1), label='p2_reordered_lpdft')
    call pmat_rasscf(p2, p2_reordered)

    do isym=1,nsym
      nao = nash(isym)
      if (nao == 0) cycle
      nio = nish(isym)
      no = norb(isym)
      nuvx = (istorp(isym+1)-istorp(isym))/nao
      do nat=1,nao
        do lrow=1,nuvx
          puvx_idx = istord(isym)+no*(lrow-1)+nio+nat
          p_idx = istorp(isym)+lrow+nuvx*(nat-1)
          energy = energy+puvx(puvx_idx)*p2_reordered(p_idx)
        end do
      end do
    end do

    energy = 0.5D0*energy
    call mma_deallocate(p2_reordered)
  end function contract_puvx_2rdm

  !***********************************************************************
  !  compute_h1e_lpdft
  !
  !> @brief Get the 1e L-PDFT Hamiltonian in the active space in MO basis.
  !> h1e = hpq + veff1 + vj + vhf_c
  !> @details
  !> @param[out] h1e_lpdft_act: 1e L-PDFT Hamiltonian in the active space in MO basis.
  !***********************************************************************
  subroutine compute_h1e_lpdft(CMO)
    real(kind=wp), dimension(ntot2), intent(in) :: CMO
    real(kind=wp), allocatable, dimension(:) :: veff1_act, vj_act, vjk_act

    call mma_allocate(veff1_act, nacpar, "VEFF1_ACT")
    call mma_allocate(vj_act, nacpar, "VJ_ACT_LPDFT")
    call mma_allocate(vjk_act, nacpar, "VJK_ACT_LPDFT")

    h1e_lpdft_ao(:) = 0.0D0
    h1e_lpdft_act(:) = 0.0D0
    veff1_act(:) = 0.0D0
    vj_act(:) = 0.0D0
    vjk_act(:) = 0.0D0

    ! Compute the hcore term and transform it to the active space in MO basis
    call get_hcore(h1e_lpdft_ao)
    call ao2moact_1e(CMO, h1e_lpdft_ao, h1e_lpdft_act, nSym, nBas, nAsh, nFro, nIsh)

    ! Add the classical Coulomb and veff2 core-response terms.
    ! Term-2: Add the classical Coulomb term
    call ao2moact_1e(CMO, vj, vj_act, nSym, nBas, nAsh, nFro, nIsh)
    call daxpy_(nacpar, 1.0D0, vj_act, 1, h1e_lpdft_act, 1)

    ! Term-4: Add the FI_V/vhf_c term, already stored in the MO basis.
    call mo2act_1e(vjk, vjk_act, nSym, nOrb, nAsh, nIsh)
    call daxpy_(nacpar, 1.0D0, vjk_act, 1, h1e_lpdft_act, 1)

    ! Compute the 1e on-top potential term and transform it to the active space
    ! in MO basis
    call mo2act_1e(veff1, veff1_act, nSym, nOrb, nAsh, nIsh)
    call daxpy_(nacpar, 1.0D0, veff1_act, 1, h1e_lpdft_act, 1)

    call mma_deallocate(vjk_act)
    call mma_deallocate(vj_act)
    call mma_deallocate(veff1_act)
  end subroutine compute_h1e_lpdft

  subroutine compute_h2e_lpdft()
    h2e_lpdft_act(:) = 0.0D0
    call get_tuvx(veff2, h2e_lpdft_act)
  end subroutine compute_h2e_lpdft

  !***********************************************************************
  !> construct_lpdft_ham
  !> Construct the L-PDFT Hamiltonian in the active space in MO basis.
  !> The L-PDFT Hamiltonian is given by:
  !> H_LPDFT = hconst + h1e + h2e
  !> @details
  !> For diagonal elements (iRoot == jRoot), we are using the 1e and 2e RDM.
  !> For off-diagonal elements (iRoot /= jRoot), we are using the 1e and 2e
  !> transition density matrices.
  !> @param[out] heff: L-PDFT Hamiltonian in the active space in MO basis.
  !***********************************************************************
  subroutine construct_lpdft_ham()

    use Constants, only: Zero
    use mspdft, only: heff
#ifdef _HDF5_
    use mh5, only: mh5_close_file, mh5_exists_dset, mh5_fetch_dset, mh5_is_hdf5, mh5_open_file_r
#endif

    implicit none
#include "warnings.h"
    integer(kind=iwp) :: iPair, iRoot, jRoot, nPairs, nSpinTdm2
    real(kind=wp), allocatable :: tdm1(:,:,:), tdm2(:,:), tdm2Spin(:,:)
#ifdef _HDF5_
    character(len=1024) :: currDir, fileName, project
    integer(kind=iwp) :: mh5id
    integer(kind=iwp), allocatable :: rootPairs(:,:)
#endif

    do iRoot=1,lRoots
      heff(iRoot,iRoot) = heff(iRoot,iRoot)+hconst
    end do

    nPairs = (lRoots*(lRoots+1))/2
    nSpinTdm2 = (nac*nac*(nac*nac+1))/2

#ifndef _HDF5_
    call WarningMessage(2,'LPDFT MS-PDFT requested, but HDF5 support is disabled.')
    call Quit_OnUserError()
#else
    call get_environment_variable('Project',project)
    if (len_trim(project) == 0) project = 'Noname'
    call get_environment_variable('CurrDir',currDir)
    if (len_trim(currDir) > 0) then
      fileName = trim(currDir)//'/'//trim(project)//'_tdm_4_lpdft.h5'
    else
      fileName = trim(project)//'_tdm_4_lpdft.h5'
    end if

    if (.not. mh5_is_hdf5(trim(fileName))) then
      call WarningMessage(2,'LPDFT transition density HDF5 file not found: '//trim(fileName))
      call Quit_OnUserError()
    end if

    mh5id = mh5_open_file_r(trim(fileName))
    if (.not. mh5_exists_dset(mh5id,'ACTIVE_1TDM')) then
      call WarningMessage(2,'LPDFT transition density HDF5 file does not contain ACTIVE_1TDM.')
      call Quit_OnUserError()
    end if
    if (.not. mh5_exists_dset(mh5id,'ACTIVE_2TDM')) then
      call WarningMessage(2,'LPDFT transition density HDF5 file does not contain ACTIVE_2TDM.')
      call Quit_OnUserError()
    end if
    if ((lRoots > 1) .and. (.not. mh5_exists_dset(mh5id,'ACTIVE_2TDM_SPINSUM'))) then
      call WarningMessage(2,'LPDFT transition density HDF5 file does not contain ACTIVE_2TDM_SPINSUM.')
      call Quit_OnUserError()
    end if

    call mma_allocate(tdm1,nac,nac,nPairs,Label='LPDFT_TDM1')
    call mma_allocate(tdm2,nacpr2,nPairs,Label='LPDFT_TDM2')
    if (lRoots > 1) call mma_allocate(tdm2Spin,nSpinTdm2,nPairs,Label='LPDFT_TDM2_SPINSUM')
    call mma_allocate(rootPairs,2,nPairs,Label='LPDFT_ROOT_PAIRS')

    tdm1(:,:,:) = Zero
    tdm2(:,:) = Zero
    call mh5_fetch_dset(mh5id,'ACTIVE_1TDM',tdm1)
    call mh5_fetch_dset(mh5id,'ACTIVE_2TDM',tdm2)
    if (lRoots > 1) then
      tdm2Spin(:,:) = Zero
      call mh5_fetch_dset(mh5id,'ACTIVE_2TDM_SPINSUM',tdm2Spin)
    end if

    if (mh5_exists_dset(mh5id,'ROOT_PAIRS')) then
      call mh5_fetch_dset(mh5id,'ROOT_PAIRS',rootPairs)
    else
      call default_root_pairs(rootPairs,nPairs)
    end if

    call mh5_close_file(mh5id)

    do iPair=1,nPairs
      iRoot = rootPairs(1,iPair)
      jRoot = rootPairs(2,iPair)
      if ((iRoot < 1) .or. (iRoot > lRoots) .or. (jRoot < 1) .or. (jRoot > lRoots)) then
        call WarningMessage(2,'LPDF transition density HDF5 file contains an invalid ROOT_PAIRS entry.')
        call Quit_OnUserError()
      end if
      if (iRoot == jRoot) then
        heff(iRoot,jRoot) = heff(iRoot,jRoot)+contract_lpdft_rdms(tdm1(:,:,iPair),tdm2(:,iPair))
      else
        heff(iRoot,jRoot) = heff(iRoot,jRoot)+contract_lpdft_tdms(tdm1(:,:,iPair),tdm2Spin(:,iPair))
      end if
      if (iRoot /= jRoot) heff(jRoot,iRoot) = heff(iRoot,jRoot)
    end do

    call print_lpdft_initial_info(fileName,nPairs)

    call mma_deallocate(rootPairs)
    if (lRoots > 1) call mma_deallocate(tdm2Spin)
    call mma_deallocate(tdm2)
    call mma_deallocate(tdm1)
#endif

  end subroutine construct_lpdft_ham

  subroutine print_lpdft_initial_info(fileName,nPairs)

    use mcpdft_output, only: iprglb
    use PrintLevel, only: USUAL

    implicit none
    character(len=*), intent(in) :: fileName
    integer(kind=iwp), intent(in) :: nPairs
    character(len=120) :: Line

    if (iprglb < USUAL) return

    write(u6,*)
    write(Line,'(6X,A)') 'LPDFT INITIAL INFORMATION'
    call CollapseOutput(1,Line)
    write(u6,'(6X,A)') repeat('-',len_trim(Line)-3)
    write(u6,*)
    write(u6,'(6X,A,T50,I6)') 'Number of roots',lRoots
    write(u6,'(6X,A,T50,I6)') 'Number of root-pair densities',nPairs
    write(u6,'(6X,A,T50,A)') 'Transition density file',trim(fileName)
    call CollapseOutput(0,'LPDFT INITIAL INFORMATION')
    write(u6,*)

  end subroutine print_lpdft_initial_info

  function contract_lpdft_rdms(rdm1,rdm2) result(ediagval)

    use Constants, only: Zero

    implicit none
    real(kind=wp), intent(in) :: rdm1(nac,nac), rdm2(nacpr2)
    real(kind=wp) :: ediagval
    integer(kind=iwp) :: iAct, iPacked, jAct

    ediagval = Zero
    iPacked = 0
    do iAct=1,nac
      do jAct=1,iAct
        iPacked = iPacked+1
        if (iAct == jAct) then
          ediagval = ediagval+h1e_lpdft_act(iPacked)*rdm1(iAct,jAct)
        else
          ediagval = ediagval+h1e_lpdft_act(iPacked)*(rdm1(iAct,jAct)+rdm1(jAct,iAct))
        end if
      end do
    end do

    ediagval = ediagval + contract_puvx_2rdm(veff2,rdm2)

  end function contract_lpdft_rdms

  function contract_spin_sum_2tdm(h2,tdm2) result(energy)

    implicit none
    real(kind=wp), intent(in) :: h2(nacpr2), tdm2(:)
    real(kind=wp) :: energy
    integer(kind=iwp) :: h2Index, iAct, ij, ijPair, ijkl, jAct, kAct, kl, klPair, lAct

    energy = 0.0_wp
    do jAct=1,nac
      do iAct=1,nac
        ij = iAct+nac*(jAct-1)
        ijPair = orbital_pair_index(iAct,jAct)
        do lAct=1,nac
          do kAct=1,nac
            kl = kAct+nac*(lAct-1)
            if (kl > ij) cycle
            klPair = orbital_pair_index(kAct,lAct)
            ijkl = pair_pair_index(ij,kl)
            h2Index = pair_pair_index(ijPair,klPair)
            energy = energy+h2(h2Index)*tdm2(ijkl)
          end do
        end do
      end do
    end do
  end function contract_spin_sum_2tdm

  function orbital_pair_index(iAct,jAct) result(idx)
    implicit none
    integer(kind=iwp), intent(in) :: iAct, jAct
    integer(kind=iwp) :: idx

    if (iAct >= jAct) then
      idx = (iAct*(iAct-1))/2+jAct
    else
      idx = (jAct*(jAct-1))/2+iAct
    end if
  end function orbital_pair_index

  function pair_pair_index(ijPair,klPair) result(idx)
    implicit none
    integer(kind=iwp), intent(in) :: ijPair, klPair
    integer(kind=iwp) :: idx

    if (ijPair >= klPair) then
      idx = (ijPair*(ijPair-1))/2+klPair
    else
      idx = (klPair*(klPair-1))/2+ijPair
    end if
  end function pair_pair_index

  !***********************************************************************
  !> contract_lpdft_tdms
  !> @brief Contract the 1e and 2e L-PDFT Hamilton
  !> with the transition density matrices to get the off-diagonal elements
  !> of the effective Hamiltonian.
  !> Note, the 2e TDM is spin-summed, and the diagonal elements are scaled by 1/2.
  !> @details
  !> @param[in] tdm1: 1e transition density matrix in the active space in MO basis.
  !> @param[in] tdm2: 2e transition density matrix in the active space in MO basis.
  !> @return eoffdiagval: The off-diagonal element of the effective Hamiltonian.
  !***********************************************************************
  function contract_lpdft_tdms(tdm1,tdm2) result(eoffdiagval)

    use Constants, only: Zero

    implicit none
    real(kind=wp), intent(in) :: tdm1(nac,nac), tdm2(:)
    real(kind=wp) :: eoffdiagval, value_1e, value_2e
    integer(kind=iwp) :: iAct, iPacked, jAct

    eoffdiagval = Zero
    value_1e = Zero
    iPacked = 0
    do iAct=1,nac
      do jAct=1,iAct
        iPacked = iPacked+1
        if (iAct == jAct) then
          value_1e = value_1e+h1e_lpdft_act(iPacked)*tdm1(iAct,jAct)
        else
          value_1e = value_1e+h1e_lpdft_act(iPacked)*(tdm1(iAct,jAct)+tdm1(jAct,iAct))
        end if
      end do
    end do

    value_2e = contract_spin_sum_2tdm(h2e_lpdft_act,tdm2)
    eoffdiagval = -(value_1e+value_2e)

  end function contract_lpdft_tdms

  subroutine default_root_pairs(rootPairs,nPairs)

    implicit none
    integer(kind=iwp), intent(in) :: nPairs
    integer(kind=iwp), intent(out) :: rootPairs(2,nPairs)
    integer(kind=iwp) :: iPair, iRoot, jRoot

    iPair = 0
    do iRoot=1,lRoots
      do jRoot=1,iRoot
        iPair = iPair+1
        rootPairs(1,iPair) = iRoot
        rootPairs(2,iPair) = jRoot
      end do
    end do

  end subroutine default_root_pairs

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
    call mma_allocate(h1e_lpdft_act, nacpar, "H1E_LPDFT_ACT")
    call mma_allocate(h2e_lpdft_act, nacpr2, "H2E_LPDFT_ACT")
    call mma_allocate(vj, ntot1, "VJ_LPDFT")
    call mma_allocate(vjk, ntot1, "VJK_LPDFT")

    call zero_global_arrays()

    dm1_core(:) = 0.0D0
    call Get_D1I_RASSCF(CMO, dm1_core)

  end subroutine init_vars

  subroutine zero_global_arrays()
    hconst = 0.0D0
    e_ot_lpdft = 0.0D0
    e_nuc_lpdft = 0.0D0
    e_wfn_lpdft = 0.0D0
    casd1_0(:) = 0.0D0
    casd1s_0(:) = 0.0D0
    casd2_0(:) = 0.0D0
    dm1_core(:) = 0.0D0
    veff1(:) = 0.0D0
    veff2(:) = 0.0D0
    h1e_lpdft_ao(:) = 0.0D0
    h1e_lpdft_act(:) = 0.0D0
    h2e_lpdft_act(:) = 0.0D0
    vj(:) = 0.0D0
    vjk(:) = 0.0D0
  end subroutine zero_global_arrays

  subroutine release_vars()
    call mma_deallocate(casd1_0)
    call mma_deallocate(casd1s_0)
    call mma_deallocate(casd2_0)
    call mma_deallocate(h1e_lpdft_ao)
    call mma_deallocate(h1e_lpdft_act)
    call mma_deallocate(h2e_lpdft_act)
    call mma_deallocate(dm1_core)
    call mma_deallocate(veff1)
    call mma_deallocate(veff2)
    call mma_deallocate(vj)
    call mma_deallocate(vjk)
  end subroutine release_vars

end module lpdft
