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
      subroutine cre_raswfn
*     SVC: Create a wavefunction file. If another .wfn file already
*     exists, it will be overwritten.
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      implicit none
#ifdef _HDF5_
#  include "rasdim.fh"
#  include "rasscf.fh"
#  include "WrkSpc.fh"
#  include "general.fh"
#  include "stdalloc.fh"
#  include "raswfn.fh"
#  include "gugx.fh"
#  include "gas.fh"
#  include "input_ras.fh"
      Integer         IDXCI(mxAct), IDXSX(mxAct)
      Common /IDSXCI/ IDXCI,        IDXSX

      integer :: dsetid
      integer, dimension(mxsym) :: NTMP1, NTMP2, NTMP3
      character(1), allocatable :: typestring(:)

*     create a new wavefunction file!
      wfn_fileid = mh5_create_file('RASWFN')

*     set module type
      call mh5_init_attr (wfn_fileid,'MOLCAS_MODULE', 'RASSCF')

*     copy basic molecular information to the HDF5 file
      call run2h5_molinfo(wfn_fileid)
      call one2h5_ovlmat(wfn_fileid, nsym, nbas)
      call one2h5_fckint(wfn_fileid, nsym, nbas)
      call one2h5_crtmom(wfn_fileid, nsym, nbas)

*     set wavefunction type
      if (iDoGAS) then
        call mh5_init_attr (wfn_fileid,'CI_TYPE', 'GAS')
      else if (IFCAS.EQ.0) then
        call mh5_init_attr (wfn_fileid,'CI_TYPE', 'CAS')
      else
        call mh5_init_attr (wfn_fileid,'CI_TYPE', 'RAS')
      end if

*     general wavefunction attributes
      call mh5_init_attr (wfn_fileid,'SPINMULT', iSpin)
      call mh5_init_attr (wfn_fileid,'LSYM', lSym)
      call mh5_init_attr (wfn_fileid,'NACTEL', nActEl)
      call mh5_init_attr (wfn_fileid,'NHOLE1', nHole1)
      call mh5_init_attr (wfn_fileid,'NELEC3', nElec3)
      call mh5_init_attr (wfn_fileid,'NCONF',  nConf)
      call mh5_init_attr (wfn_fileid,'NSTATES', nRoots)
      call mh5_init_attr (wfn_fileid,'NROOTS', lRoots)

      call mh5_init_attr (wfn_fileid,'L2ACT', 1, [mxAct], IDXSX)
      call mh5_init_attr (wfn_fileid,'A2LEV', 1, [mxAct], IDXCI)

*     iteration(s)
      wfn_iter = mh5_create_attr_int (wfn_fileid,'RASSCF_ITERATIONS')

*     molecular orbital type index
      if (iDoGAS) then
        NTMP1(:)=0
        NTMP2(:)=sum(NGSSH(1:NGAS,:),dim=1)
        NTMP3(:)=0
      else
        NTMP1(:)=NRS1(:)
        NTMP2(:)=NRS2(:)
        NTMP3(:)=NRS3(:)
      end if
      call mma_allocate(typestring, ntot)
      call orb2tpstr(NSYM,NBAS,
     $        NFRO,NISH,NTMP1,NTMP2,NTMP3,NSSH,NDEL,
     $        typestring)
      dsetid = mh5_create_dset_str(wfn_fileid,
     $        'MO_TYPEINDICES', 1, [NTOT],1)
      call mh5_init_attr(dsetid, 'description',
     $        'Type index of the molecular orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')
      call mh5_put_dset(dsetid, typestring)
      call mma_deallocate(typestring)
      call mh5_close_dset(dsetid)

*     roots
      call mh5_init_attr (wfn_fileid,
     $        'STATE_ROOTID', 1, [nRoots], iRoot)
      call mh5_init_attr (wfn_fileid,
     $        'STATE_WEIGHT', 1, [nRoots], Weight)

*     energy (for each CI root)
      wfn_energy = mh5_create_dset_real (wfn_fileid,
     $        'ROOT_ENERGIES', 1, [lRoots])
      call mh5_init_attr(wfn_energy, 'description',
     $        'Energy for each root in the CI, '//
     $        'arranged as array of [NROOTS]')

*     molecular orbital coefficients
      wfn_mocoef = mh5_create_dset_real(wfn_fileid,
     $        'MO_VECTORS', 1, [NTOT2])
      call mh5_init_attr(wfn_mocoef, 'description',
     $        'Coefficients of the average orbitals, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

*     molecular orbital occupation numbers
      wfn_occnum = mh5_create_dset_real(wfn_fileid,
     $        'MO_OCCUPATIONS', 1, [NTOT])
      call mh5_init_attr(wfn_occnum, 'description',
     $        'Occupation numbers of the average orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')

*     molecular orbital energies
      wfn_orbene = mh5_create_dset_real(wfn_fileid,
     $        'MO_ENERGIES', 1, [NTOT])
      call mh5_init_attr(wfn_orbene, 'description',
     $        'Orbital energies of the average orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')

*     molecular orbital symmetry irreps
      wfn_supsym = mh5_create_dset_int(wfn_fileid,
     $        'SUPSYM_IRREP_INDICES', 1, [NTOT])
      call mh5_init_attr(wfn_supsym, 'description',
     $        'Supersymmetry ID of the average orbitals '//
     $        'arranged as blocks of size [NBAS(i)], i=1,#irreps')
      call mh5_init_attr(wfn_supsym, 'pointgroup', 'D5h')
      call mh5_init_attr(wfn_supsym, '#irreps', 7)
      call mh5_put_dset(wfn_supsym, IXSYM)

*     CI data for each root
      wfn_cicoef = mh5_create_dset_real(wfn_fileid,
     $        'CI_VECTORS', 2, [nConf, lRoots])
      call mh5_init_attr(wfn_cicoef, 'description',
     $        'Coefficients of configuration state functions '//
     $        'in Split-GUGA ordering, size [NCONF] '//
     $        'for each root in NROOTS: [NROOTS,NCONF].')

*     density matrices for each root
      wfn_dens = mh5_create_dset_real (wfn_fileid,
     $        'DENSITY_MATRIX', 3, [NAC, NAC, lRoots])
      call mh5_init_attr(wfn_dens, 'description',
     $        'active 1-body density matrix, size [NAC,NAC] '//
     $        'for each root in NROOTS: [NROOTS,NAC,NAC].')

      wfn_spindens = mh5_create_dset_real (wfn_fileid,
     $        'SPINDENSITY_MATRIX', 3, [NAC, NAC, lRoots])
      call mh5_init_attr(wfn_spindens, 'description',
     $        'active 1-body spin density matrix, size [NAC,NAC] '//
     $        'for each root in NROOTS: [NROOTS,NAC,NAC].')

      if (KeyTDM) then
      wfn_transdens = mh5_create_dset_real (wfn_fileid,
     $        'TRANSITION_DENSITY_MATRIX', 3,
     $        [NAC, NAC, lRoots*(lRoots-1)/2])
      call mh5_init_attr(wfn_transdens, 'description',
     $        'active 1-body transition density matrix, '//
     $        'size [NAC,NAC] for each pair of roots in NROOTS: '//
     $        '[NROOTS*(NROOTS-1)/2,NAC,NAC].')
      end if

#ifdef _DMRG_
      if (doDMRG) then
! Leon 1/12/2016: Add the QCMaquis checkpoint name to the description of each state
! maximum allowed filename length is equal to MH5_MAX_LBL_LEN=256
        wfn_dmrg_checkpoint = mh5_create_dset_str(wfn_fileid,
     $        'QCMAQUIS_CHECKPOINT', 1, [lRoots], 256)
        call mh5_init_attr(wfn_dmrg_checkpoint,'description',
     $        'QCMaquis checkpoint directory names for each root'//
     $        ' in [NROOTS].')
      end if
#endif
#endif
      end
