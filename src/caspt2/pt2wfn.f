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
*               2018, Ignacio Fdez. Galvan                             *
************************************************************************
      module pt2wfn
      integer :: pt2wfn_id
      logical :: pt2wfn_is_h5 = .False.
      integer :: pt2wfn_refene, pt2wfn_energy
      integer :: pt2wfn_mocoef, pt2wfn_occnum, pt2wfn_orbene
      integer :: pt2wfn_cicoef
      integer :: pt2wfn_heff
      integer :: pt2wfn_dens
      save

      contains

      subroutine pt2wfn_init
*     SVC: Create a wavefunction file. If another .wfn file already
*     exists, it will be overwritten.
      use refwfn
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "pt2_guga.fh"
#ifdef _HDF5_
#  include "mh5.fh"

      integer :: dsetid, ndmat, i
      character(1), allocatable :: typestring(:)
#endif

      If (refwfn_active) Then
        Call WarningMessage(2,'Active reference wavefunction file, '//
     &    'cannot create new PT2 wavefunction file, aborting!')
        Call AbEnd
      End If

#ifdef _HDF5_
      If (refwfn_is_h5) Then
        pt2wfn_is_h5 = .True.
*     create a new wavefunction file!
        pt2wfn_id = mh5_create_file('PT2WFN')
*     no JOBMIX for HDF5 files
        ifmix = .False.

*     set module type
        call mh5_init_attr (pt2wfn_id,'MOLCAS_MODULE', 'CASPT2')

*     copy basic molecular information to the HDF5 file
        call run2h5_molinfo(pt2wfn_id)
        call one2h5_ovlmat(pt2wfn_id, nsym, nbas)
        call one2h5_fckint(pt2wfn_id, nsym, nbas)
        call one2h5_crtmom(pt2wfn_id, nsym, nbas)

*     set wavefunction type
        if (NRAS1T+NRAS3T.EQ.0) then
          call mh5_init_attr (pt2wfn_id,'CI_TYPE', 'CAS')
        else
          call mh5_init_attr (pt2wfn_id,'CI_TYPE', 'RAS')
        end if

*     general wavefunction attributes
        call mh5_init_attr (pt2wfn_id,'SPINMULT', iSpin)
        call mh5_init_attr (pt2wfn_id,'LSYM', lSym)
        call mh5_init_attr (pt2wfn_id,'NACTEL', nActEl)
        call mh5_init_attr (pt2wfn_id,'NHOLE1', nHole1)
        call mh5_init_attr (pt2wfn_id,'NELEC3', nEle3)
        call mh5_init_attr (pt2wfn_id,'NCONF',  nConf)
        call mh5_init_attr (pt2wfn_id,'NSTATES', NSTATE)

        call mh5_init_attr (pt2wfn_id,'L2ACT', 1, [mxAct], L2ACT)
        call mh5_init_attr (pt2wfn_id,'A2LEV', 1, [mxAct], LEVEL)

*     molecular orbital type index
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

*     roots
        call mh5_init_attr (pt2wfn_id,
     $        'STATE_ROOTID', 1, [NSTATE], MSTATE)
        call mh5_init_attr (pt2wfn_id,
     $        'ROOT2STATE', 1, [LROOTS], ROOT2STATE)

*     reference energy (for each CI root)
        pt2wfn_refene = mh5_create_dset_real (pt2wfn_id,
     $        'STATE_REFWF_ENERGIES', 1, [NSTATE])
        call mh5_init_attr(pt2wfn_refene, 'description',
     $        'Reference energy for each state, '//
     $        'arranged as array of [NSTATES]')

*     PT2 energy (for each CI root)
        pt2wfn_energy = mh5_create_dset_real (pt2wfn_id,
     $        'STATE_PT2_ENERGIES', 1, [NSTATE])
        call mh5_init_attr(pt2wfn_energy, 'description',
     $        'PT2 energy for each state, '//
     $        'arranged as array of [NSTATES]')

*     molecular orbital coefficients
        pt2wfn_mocoef = mh5_create_dset_real(pt2wfn_id,
     $        'MO_VECTORS', 1, [NBSQT])
        call mh5_init_attr(pt2wfn_mocoef, 'description',
     $        'Coefficients of the average orbitals, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

*     molecular orbital occupation numbers
*     (most probably empty, but left for compatibility)
        pt2wfn_occnum = mh5_create_dset_real(pt2wfn_id,
     $        'MO_OCCUPATIONS', 1, [NBAST])
        call mh5_init_attr(pt2wfn_occnum, 'description',
     $        'Occupation numbers of the average orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')

*     molecular orbital energies
*     (most probably empty, but left for compatibility)
        pt2wfn_orbene = mh5_create_dset_real(pt2wfn_id,
     $        'MO_ENERGIES', 1, [NBAST])
        call mh5_init_attr(pt2wfn_orbene, 'description',
     $        'Orbital energies of the average orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')

*     CI data for each root
        pt2wfn_cicoef = mh5_create_dset_real(pt2wfn_id,
     $        'CI_VECTORS', 2, [nConf, NSTATE])
        call mh5_init_attr(pt2wfn_cicoef, 'description',
     $        'Coefficients of configuration state functions '//
     $        'in Split-GUGA ordering for each STATE, '//
     $        'arranged as matrix of size [NCONF,NSTATES]')

*     effective Hamiltonian coefficients
        If (IFMSCOUP) Then
          pt2wfn_heff = mh5_create_dset_real(pt2wfn_id,
     $        'H_EFF', 2, [NSTATE, NSTATE])
          call mh5_init_attr(pt2wfn_heff, 'description',
     $        'Effective (X)MS-CASPT2 Hamiltonian, '//
     $        'arranged as matrix of size [NSTATES,NSTATES]')
        End If

*     density matrices
        If (IFPROP) Then
          ndmat=0
          Do i=1,NSYM
            ndmat=ndmat+(NORB(i)**2+NORB(i))/2
          End Do

          pt2wfn_dens = mh5_create_dset_real(pt2wfn_id,
     $        'DENSITY_MATRIX', 2, [ndmat, NSTATE])
          call mh5_init_attr(pt2wfn_dens, 'description',
     $        '1-body density matrix, arranged as blocks of size '//
     $        'NDMAT=sum([NORB(i)*(NORB(i)+1)/2], i=1,#irreps), '//
     $        'where NORB excludes frozen and deleted orbitals, '//
     $        'for each state: [NDMAT,NSTATES].')
        End If

      Else
#endif
        pt2wfn_is_h5 = .False.
#ifdef _HDF5_
      End If
#endif
      end subroutine

      subroutine pt2wfn_data
      use refwfn
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      real*8, allocatable :: BUF(:)
      integer :: ISTATE, IDISK

      If (pt2wfn_is_h5) Then
        call mma_allocate(BUF,NCONF)
        IDISK = IDCIEX
        DO ISTATE=1,NSTATE
          CALL DDAFILE(LUCIEX,2,BUF,NCONF,IDISK)
          call mh5_put_dset_array_real(pt2wfn_cicoef,
     &           BUF,[NCONF,1],[0,ISTATE-1])
        END DO
        call mma_deallocate(BUF)

        call mma_allocate(BUF,NCMO)
        IDISK = IAD1M(1)
        CALL DDAFILE(LUONEM,2,BUF,NCMO,IDISK)
        call mh5_put_dset_array_real(pt2wfn_mocoef,BUF)
        call mma_deallocate(BUF)
      End If
#endif
      end subroutine

      subroutine pt2wfn_estore(Heff)
      use refwfn
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
      real*8 :: Heff(nstate,nstate)
#ifdef _HDF5_
#  include "mh5.fh"
      If (pt2wfn_is_h5) Then
        call mh5_put_dset_array_real(pt2wfn_energy, ENERGY)
        call mh5_put_dset_array_real(pt2wfn_refene, REFENE)
        If (IFMSCOUP) Then
          call mh5_put_dset_array_real(pt2wfn_heff, Heff)
        End If
      End If
#else
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Heff)
#endif
      end subroutine

      subroutine pt2wfn_densstore(Dmat,nDmat)
      use refwfn
      implicit none
#include "rasdim.fh"
#include "caspt2.fh"
      integer :: nDmat
      real*8 :: Dmat(nDmat)
#ifdef _HDF5_
#  include "mh5.fh"
      If (pt2wfn_is_h5) Then
        call mh5_put_dset_array_real(pt2wfn_dens, Dmat,
     $                               [nDmat, 1], [0, JSTATE-1])
      End If
#else
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Dmat)
#endif
      end subroutine

      subroutine pt2wfn_close
#ifdef _HDF5_
      if (pt2wfn_is_h5) then
        call mh5_close_file(pt2wfn_id)
      end if
      pt2wfn_id = -1
#endif
      end subroutine
      end module
