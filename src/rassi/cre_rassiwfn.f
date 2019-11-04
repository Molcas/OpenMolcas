************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine cre_rassiwfn
*     SVC: Create a wavefunction file. If another .wfn file already
*     exists, it will be overwritten.
      use kVectors
#ifdef _HDF5_
      use rassi_global_arrays, only: JBNUM
#endif
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      implicit none
#ifdef _HDF5_
#  include "Molcas.fh"
#  include "cntrl.fh"
#  include "rassi.fh"
#  include "symmul.fh"
#  include "WrkSpc.fh"
#  include "stdalloc.fh"
#  include "rassiwfn.fh"
#  include "lebedev.fh"

      integer :: ISTATE, NSS
      integer :: nData, nIJ
      integer, allocatable :: state_irreps(:), state_mult(:)
      integer :: nbast

      nbast = sum(nbasf(1:nsym)**2)

*     create a new wavefunction file!
      wfn_fileid = mh5_create_file('RASSIWFN')

*     set module type
#ifdef _DMRG_
      if(doDMRG)then
        call mh5_init_attr (wfn_fileid,'MOLCAS_MODULE', 'MPSSI')
      else
#endif
      call mh5_init_attr (wfn_fileid,'MOLCAS_MODULE', 'RASSI')
#ifdef _DMRG_
      end if
#endif

*     copy basic molecular information to the HDF5 file
      call run2h5_molinfo(wfn_fileid)
      call one2h5_ovlmat(wfn_fileid, nsym, nbasf)
      call one2h5_crtmom(wfn_fileid, nsym, nbasf)

*     general wavefunction attributes
      call mh5_init_attr (wfn_fileid,'NSTATE', NSTATE)

      NSS=0
      DO ISTATE=1,NSTATE
        NSS=NSS+MLTPLT(JBNUM(ISTATE))
      END DO

*     irrep per state
      call mma_allocate(state_irreps, NSTATE)
      do istate=1,nstate
        state_irreps(istate) = IRREP(JBNUM(ISTATE))
      end do
      call mh5_init_attr (wfn_fileid,
     $        'STATE_IRREPS', 1, [NSTATE], state_irreps)
      call mma_deallocate(state_irreps)

*     multiplicity per state
      call mma_allocate(state_mult, NSTATE)
      do istate=1,nstate
        state_mult(istate) = MLTPLT(JBNUM(ISTATE))
      end do
      call mh5_init_attr (wfn_fileid,
     $        'STATE_SPINMULT', 1, [NSTATE], state_mult)
      call mma_deallocate(state_mult)

*     overlaps of the input states
      wfn_overlap = mh5_create_dset_real(wfn_fileid,
     $        'ORIGINAL_OVERLAPS', 2, [NSTATE,NSTATE])
      call mh5_init_attr(wfn_overlap, 'description',
     $        'Overlaps between the original (input) states, '//
     $        'a symmetric matrix of size [NSTATE,NSTATE]')

*     energies of the orthogonal (output) spin free states (SFS)
      wfn_sfs_energy = mh5_create_dset_real (wfn_fileid,
     $        'SFS_ENERGIES', 1, [NSTATE])
      call mh5_init_attr(wfn_sfs_energy, 'description',
     $        'Energy for each spin-free state, '//
     $        'arranged as array of size [NSTATE]')

*     SFS coefficients
      wfn_sfs_coef = mh5_create_dset_real(wfn_fileid,
     $        'SFS_COEFFICIENTS', 2, [NSTATE,NSTATE])
      call mh5_init_attr(wfn_sfs_coef, 'description',
     $        'Eigenstates of the input Hamiltonian, '//
     $        'expressed as linear combinations of the input '//
     $        'states, as a matrix of size [NSTATE,NSTATE], '//
     $        'fast index corresponds to input states.')

*     SFS properties
      wfn_sfs_angmom = mh5_create_dset_real(wfn_fileid,
     $        'SFS_ANGMOM', 3, [NSTATE,NSTATE,3])
      call mh5_init_attr(wfn_sfs_angmom, 'description',
     $        'Angular momentum components between the spin-free '//
     $        'states stored as <SFS1|iL(x,y,z)|SFS2> in'//
     $        ' [3,NSTATE,NSTATE]')

      wfn_sfs_edipmom = mh5_create_dset_real(wfn_fileid,
     $        'SFS_EDIPMOM', 3, [NSTATE,NSTATE,3])
      call mh5_init_attr(wfn_sfs_edipmom, 'description',
     $        'Electric dipole momentum components between the '//
     $        'spin-free states stored as <SFS1|ED(x,y,z)|SFS2> in'//
     $        ' [3,NSTATE,NSTATE]')

      wfn_sfs_amfi = mh5_create_dset_real(wfn_fileid,
     $        'SFS_AMFIINT', 3, [NSTATE,NSTATE,3])
      call mh5_init_attr(wfn_sfs_amfi, 'description',
     $        'Components of the spin-orbit integrals between the '//
     $        'spin-free states stored as '//
     $        '<SFS1|spin-orbit-operator|SFS2> in'//
     $        ' [3,NSTATE,NSTATE]')

*     SFS transition density
      wfn_sfs_tdm = mh5_create_dset_real(wfn_fileid,
     $        'SFS_TRANSITION_DENSITIES', 3, [NBAST,NSTATE,NSTATE])
      call mh5_init_attr(wfn_sfs_tdm, 'description',
     $        'Transition density matrices for each pair of spin-'//
     $        'free states, matrix of size [NSTATE,NSTATE,NBAST], '//
     $        'where NBAST is the sum of NBAS(I)**2 for I=1,NSYM.'//
     $        'Only contributing symmetry blocks are stored')

*     SFS spin transition density
      wfn_sfs_tsdm = mh5_create_dset_real(wfn_fileid,
     $        'SFS_TRANSITION_SPIN_DENSITIES', 3, [NBAST,NSTATE,NSTATE])
      call mh5_init_attr(wfn_sfs_tsdm, 'description',
     $        'Transition spin density matrices for each pair of spin'//
     $        '-free states, matrix of size [NSTATE,NSTATE,NBAST], '//
     $        'where NBAST is the sum of NBAS(I)**2 for I=1,NSYM.'//
     $        'Only contributing symmetry blocks are stored')

      if (do_tmom) then
*     SFS intermediate transition vectors
      nIJ=NSTATE*(NSTATE-1)/2
      nData= 1 + 3 + 2*3
      wfn_sfs_tm = mh5_create_dset_real(wfn_fileid,
     $        'SFS_TRANSITION_VECTORS', 4, [nData,nQuad,nIJ,nk_Vector])
      call mh5_init_attr(wfn_sfs_tm, 'description',
     $        'SFS intermediate transition vectors '//
     $        'for nVec specific k-vectors, '//
     $        'nIJ unique pairs of SF states, '//
     $        'and nQuad quadrature k-vectors; '//
     $        'stored as weight, k-vector (3 components), '//
     $        'real transition vector (3 components), '//
     $        'and imaginary transition vector (3 components), '//
     $        'giving a matrix of size [nVec,nIJ,nQuad,nData]')
      end if

      if (ifso) then

*     SFS WE-reduced TDMs of triplet type
      wfn_sfs_wetdm = mh5_create_dset_real(wfn_fileid,
     $        'SFS_WE_TRANSITION_DENSITIES', 3, [NBAST,NSTATE,NSTATE])
      call mh5_init_attr(wfn_sfs_wetdm, 'description',
     $        'WE-reduced TDMs for each pair of spin-free states,'//
     $        'matrix of size [NSTATE,NSTATE,NBAST], where NBAST '//
     $        'is the sum of NBAS(I)**2 for I=1,NSYM.'//
     $        'Only contributing symmetry blocks are stored')

*     energies of the spin orbit states (SOS)
      wfn_sos_energy = mh5_create_dset_real (wfn_fileid,
     $        'SOS_ENERGIES', 1, [NSS])
      call mh5_init_attr(wfn_sos_energy, 'description',
     $        'Energy for each spin-orbit state, '//
     $        'arranged as array of [NSS]')

*     SO complex Hamiltonian
      wfn_sos_hsor = mh5_create_dset_real(wfn_fileid,
     $        'HSO_MATRIX_REAL', 2, [NSS,NSS])
      call mh5_init_attr(wfn_sos_hsor, 'description',
     $        'The spin-orbit Hamiltonian, '//
     $        '2D-array, real part as [NSS,NSS]')
      wfn_sos_hsoi = mh5_create_dset_real(wfn_fileid,
     $        'HSO_MATRIX_IMAG', 2, [NSS,NSS])
      call mh5_init_attr(wfn_sos_hsoi, 'description',
     $        'The spin-orbit Hamiltonian, '//
     $        '2D-array, imaginary part as [NSS,NSS]')

*     SOS coefficients
      wfn_sos_coefr = mh5_create_dset_real(wfn_fileid,
     $        'SOS_COEFFICIENTS_REAL', 2, [NSS,NSS])
      call mh5_init_attr(wfn_sos_coefr, 'description',
     $        'Eigenstates of the spin-orbit Hamiltonian, '//
     $        'expressed as linear combinations of the spin-free '//
     $        'states, 2D-array of real part as [NSS,NSS], '//
     $        'fast index corresponds to spin-free states.')
      wfn_sos_coefi = mh5_create_dset_real(wfn_fileid,
     $        'SOS_COEFFICIENTS_IMAG', 2, [NSS,NSS])
      call mh5_init_attr(wfn_sos_coefi, 'description',
     $        'Eigenstates of the spin-orbit Hamiltonian, '//
     $        'expressed as linear combinations of the spin-free '//
     $        'states, 2D-array of imaginary part as [NSS,NSS], '//
     $        'fast index corresponds to spin-free states.')

*     SOS properties
      wfn_sos_angmomr = mh5_create_dset_real(wfn_fileid,
     $        'SOS_ANGMOM_REAL', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_angmomr, 'description',
     $        'Angular momentum components between the spin-orbit '//
     $        'states stored as <SOS1|iL(x,y,z)|SOS2> in'//
     $        ' [3,NSS,NSS], real part')

      wfn_sos_angmomi = mh5_create_dset_real(wfn_fileid,
     $        'SOS_ANGMOM_IMAG', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_angmomi, 'description',
     $        'Angular momentum components between the spin-orbit '//
     $        'states stored as <SOS1|iL(x,y,z)|SOS2> in'//
     $        ' [3,NSS,NSS], imaginary part')

      wfn_sos_spinr = mh5_create_dset_real(wfn_fileid,
     $        'SOS_SPIN_REAL', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_spinr, 'description',
     $        'Spin momentum components between the spin-orbit '//
     $        'states stored as <SOS1|SPIN(x,y,z)|SOS2> in'//
     $        ' [NSS,NSS,3], real part')

      wfn_sos_spini = mh5_create_dset_real(wfn_fileid,
     $        'SOS_SPIN_IMAG', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_spini, 'description',
     $        'Spin momentum components between the spin-orbit '//
     $        'states stored as <SOS1|SPIN(x,y,z)|SOS2> in'//
     $        ' [NSS,NSS,3], imaginary part')

      wfn_sos_edipmomr = mh5_create_dset_real(wfn_fileid,
     $        'SOS_EDIPMOM_REAL', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_edipmomr, 'description',
     $        'Electric dipole momentum components between the '//
     $        'spin-orbit states stored as <SOS1|ED(x,y,z)|SOS2> in'//
     $        ' [3,NSS,NSS], real part')

      wfn_sos_edipmomi = mh5_create_dset_real(wfn_fileid,
     $        'SOS_EDIPMOM_IMAG', 3, [NSS,NSS,3])
      call mh5_init_attr(wfn_sos_edipmomi, 'description',
     $        'Electric dipole momentum components between the '//
     $        'spin-orbit states stored as <SOS1|ED(x,y,z)|SOS2> in'//
     $        ' [3,NSS,NSS], imaginary part')

      if (do_tmom) then
      nIJ=NSS*(NSS-1)/2
      nData= 1 + 3 + 2*3
      wfn_sos_tm = mh5_create_dset_real(wfn_fileid,
     $        'SOS_TRANSITION_VECTORS', 4,
     $        [nData,2*nQuad,nIJ,nk_Vector])
      call mh5_init_attr(wfn_sos_tm, 'description',
     $        'SOS intermediate transition vectors '//
     $        'for nVec specific k-vectors, '//
     $        'nIJ unique pairs of SF states, '//
     $        'and nQuad quadrature k-vectors; '//
     $        'stored as weight, k-vector (3 components), '//
     $        'real transition vector (3 components), '//
     $        'and imaginary transition vector (3 components), '//
     $        'giving a matrix of size [nVec,nIJ,nQuad,nData]')
      end if

      end if
#endif
      end
