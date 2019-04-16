************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016, Steven Vancoillie                                *
*               2017, Stefan Knecht                                    *
************************************************************************
! author: S. Knecht (modified version of pt2wfn.f in CASPT2 written by
!                    S. Vancoillie)
!***********************************************************************
      module nevpt2wfn
      integer :: pt2wfn_id
      logical :: pt2wfn_is_h5 = .False.
      integer :: pt2wfn_refene, pt2wfn_energy_sc, pt2wfn_energy_pc
      integer :: pt2wfn_mocoef!, pt2wfn_occnum, pt2wfn_orbene
      integer :: pt2wfn_heff_sc, pt2wfn_heff_pc
      integer :: pt2wfn_heff_evc_sc, pt2wfn_heff_evc_pc
      integer :: pt2wfn_ref_checkpoint
      save

      contains

      subroutine nevpt2wfn_init(create_h5)
!     Create a wavefunction file and replace any existing .wfn file
      use refwfn
      use nevpt2_cfg
      use info_state_energy  ! energies, effective Hamiltonian
      use info_orbital_space ! orbital specifications read from JobIph
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "pt2_guga.fh"
      logical, intent(in)       :: create_h5
#ifdef _HDF5_
#  include "mh5.fh"
      integer :: dsetid, i
      character(1), allocatable :: typestring(:)
#endif

      If (refwfn_active) Then
        Call WarningMessage(2,'Active reference wavefunction file, '//
     &    'cannot create new PT2 wavefunction file, aborting!')
        Call AbEnd
      End If

#ifdef _HDF5_
      If (create_h5) Then
        pt2wfn_is_h5 = .True.

      !> create a new wavefunction file!
        pt2wfn_id = mh5_create_file('NEVPT2WFN')

      !> set module type
        call mh5_init_attr (pt2wfn_id,'MOLCAS_MODULE', 'NEVPT2')

      !> copy basic molecular information to the HDF5 file
        call run2h5_molinfo(pt2wfn_id)
        call one2h5_ovlmat(pt2wfn_id, nsym, nbas)
        call one2h5_fckint(pt2wfn_id, nsym, nbas)
        call one2h5_crtmom(pt2wfn_id, nsym, nbas)

      !> set wavefunction type
        call mh5_init_attr (pt2wfn_id,'CI_TYPE', 'CAS')

      !> general wavefunction attributes
        call mh5_init_attr (pt2wfn_id,'SPINMULT', nSpin)
        call mh5_init_attr (pt2wfn_id,'LSYM', lSym)
        call mh5_init_attr (pt2wfn_id,'NACTEL', nr_active_electrons)
        call mh5_init_attr (pt2wfn_id,'NHOLE1', 0)
        call mh5_init_attr (pt2wfn_id,'NELEC3', 0)
        call mh5_init_attr (pt2wfn_id,'NCONF',  1)
        call mh5_init_attr (pt2wfn_id,'NSTATES', nr_states)

      !> keep it for compatibility
        call mh5_init_attr (pt2wfn_id,'L2ACT', 1, [MXLEV], L2ACT)
        call mh5_init_attr (pt2wfn_id,'A2LEV', 1, [MXLEV], LEVEL)

      !> molecular orbital type index
        call mma_allocate(typestring, nbast)
        call orb2tpstr(NSYM,NBAS,
     $        NFRO,NISH,NRAS1,NRAS2,NRAS3,NSSH,NDEL,
     $        typestring)
        dsetid = mh5_create_dset_str(pt2wfn_id,
     $        'MO_TYPEINDICES', 1, [NBAST],1)
        call mh5_init_attr(dsetid, 'description',
     $        'Type index of the molecular orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')
        call mh5_put_dset(dsetid, typestring)
        call mma_deallocate(typestring)
        call mh5_close_dset(dsetid)

      !> roots
!       call mh5_init_attr (pt2wfn_id,
!    $        'STATE_ROOTID', 1, [nr_states], MultGroup%state)

      !> setup dummy root to state translation
        root2state = 0
        do i=1,nr_states
          root2state(i)=i
        end do

        call mh5_init_attr (pt2wfn_id,
     $        'STATE_ROOTID', 1, [nr_states], root2state)

      !> setup dummy root to state translation - part 2 -
        root2state = 0
        do i=1,nr_states
          root2state(i) = MultGroup%state(i)
        end do

        call mh5_init_attr (pt2wfn_id,
     $        'ROOT2STATE', 1, [nr_states], ROOT2STATE)

      !> reference energy (for each state)
        pt2wfn_refene = mh5_create_dset_real (pt2wfn_id,
     $        'STATE_REFWF_ENERGIES', 1, [nr_states])
        call mh5_init_attr(pt2wfn_refene, 'description',
     $        'Reference energy for each state, '//
     $        'arranged as array of [nr_states]')

      !> PT2 energy (SC for each state)
      pt2wfn_energy_sc = mh5_create_dset_real (pt2wfn_id,
     $      'STATE_PT2_ENERGIES_SC', 1, [nr_states])
      call mh5_init_attr(pt2wfn_energy_sc, 'description',
     $      'PT2 energy (SC) for each state, '//
     $      'arranged as array of [nr_states]')
      !> effective Hamiltonian (SC)
        pt2wfn_heff_sc = mh5_create_dset_real(pt2wfn_id,
     &        'H_EFF_SC', 2, [nr_states, nr_states])
        call mh5_init_attr(pt2wfn_heff_sc, 'description',
     &        'Effective QD-NEVPT2 hamiltonian (SC), '//
     &        'arranged as matrix of size [nr_states,nr_states]')

      !> molecular orbital coefficients
        pt2wfn_mocoef = mh5_create_dset_real(pt2wfn_id,
     $        'MO_VECTORS', 1, [NBSQT])
        call mh5_init_attr(pt2wfn_mocoef, 'description',
     $        'Coefficients of the average orbitals, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

      if(.not.no_pc)then
        !> PT2 energy (PC for each state) - default
        pt2wfn_energy_pc = mh5_create_dset_real (pt2wfn_id,
     $        'STATE_PT2_ENERGIES', 1, [nr_states])
        call mh5_init_attr(pt2wfn_energy_pc, 'description',
     $        'PT2 energy (PC) for each state, '//
     $        'arranged as array of [nr_states]')
        !> effective Hamiltonian (PC) - default
        pt2wfn_heff_pc = mh5_create_dset_real(pt2wfn_id,
     &        'H_EFF', 2, [nr_states, nr_states])
        call mh5_init_attr(pt2wfn_heff_pc, 'description',
     &        'Effective QD-NEVPT2 hamiltonian (PC), '//
     &        'arranged as matrix of size [nr_states,nr_states]')

      end if

#ifdef _DMRG_
        !> maximum allowed filename length is equal to MH5_MAX_LBL_LEN=256
        pt2wfn_ref_checkpoint = mh5_create_dset_str(pt2wfn_id,
     $        'QCMAQUIS_CHECKPOINT', 1, [nr_states], 256)
        call mh5_init_attr(pt2wfn_ref_checkpoint,'description',
     $        'QCMaquis checkpoint directory names for each root'//
     $        ' in [nr_states].')
#endif
      Else
#endif
        pt2wfn_is_h5 = .False.
#ifdef _HDF5_
      End If
#endif
      end subroutine

      subroutine nevpt2wfn_data
#ifdef _DMRG_
      use qcmaquis_info
#endif
      use refwfn
      use nevpt2_cfg, only : MultGroup
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      real*8, allocatable :: BUF(:)
      integer :: IDISK
#ifdef _DMRG_
      integer :: i
#endif

      If (pt2wfn_is_h5) Then
        call mma_allocate(BUF,NCMO)
        IDISK = 0
        CALL DDAFILE(LUONEM,2,BUF,NCMO,IDISK)
        call mh5_put_dset_array_real(pt2wfn_mocoef,BUF)
        call mma_deallocate(BUF)
#ifdef _DMRG_
        if(allocated(qcm_group_names))then
          do i = 1, size(MultGroup%State)
            call mh5_put_dset_array_str(
     &           pt2wfn_ref_checkpoint,
     &           qcm_group_names(1)%states(MultGroup%State(i)),
     &           [1],
     &           [i-1]
     &          )
          end do
        end if
#endif

      End If
#endif
      end subroutine

      subroutine nevpt2wfn_estore()
      use refwfn
      use nevpt2_cfg
      use info_state_energy  ! energies + effective Hamiltonian
      implicit none
#ifdef _HDF5_
#  include "mh5.fh"

      If (pt2wfn_is_h5) Then
        !> reference energies aka DMRG-SCF energies
        call mh5_put_dset_array_real(pt2wfn_refene, e)
        call mh5_put_dset_array_real(pt2wfn_energy_sc, psimp)
        !> effective Hamiltonian
        call mh5_put_dset_array_real(pt2wfn_heff_sc,e2mp)
        if(.not.no_pc)then
          !> single-state PC energies
          call mh5_put_dset_array_real(pt2wfn_energy_pc, psien)
          !> effective PC Hamiltonian
          call mh5_put_dset_array_real(pt2wfn_heff_pc,e2en)
        end if
      End If
#endif
      end subroutine nevpt2wfn_estore

      subroutine nevpt2wfn_close
#ifdef _DMRG_
      use qcmaquis_info
#endif
      implicit none
#ifdef _HDF5_
      if (pt2wfn_is_h5) then
        call mh5_close_file(pt2wfn_id)
        pt2wfn_is_h5 = .False.
        pt2wfn_id    = -1
#ifdef _DMRG_
        call qcmaquis_info_deinit
#endif

      end if
#endif
      end subroutine

      end module nevpt2wfn
