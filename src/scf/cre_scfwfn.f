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
* Copyright (C) 2015,2016, Steven Vancoillie                           *
*               2018, Ignacio Fdez. Galvan                             *
************************************************************************
      subroutine cre_scfwfn
*     SVC: Create a wavefunction file. If another .scf.h5 file already
*     exists, it will be overwritten.
#ifdef _HDF5_
      use mh5, only: mh5_create_file, mh5_init_attr,
     &               mh5_create_dset_real, mh5_create_dset_str
#endif
      implicit none
#ifdef _HDF5_
#  include "mxdm.fh"
#  include "stdalloc.fh"
#  include "scfwfn.fh"
#  include "infscf.fh"

*     create a new wavefunction file!
      wfn_fileid = mh5_create_file('SCFWFN')

*     set module type
      call mh5_init_attr (wfn_fileid,'MOLCAS_MODULE', 'SCF')

*     copy basic molecular information to the HDF5 file
      call run2h5_molinfo(wfn_fileid)
      call one2h5_ovlmat(wfn_fileid, nsym, nbas)
      call one2h5_fckint(wfn_fileid, nsym, nbas)
      call one2h5_crtmom(wfn_fileid, nsym, nbas)

*     energy
      wfn_energy = mh5_create_dset_real (wfn_fileid,'ENERGY')
      call mh5_init_attr(wfn_energy, 'DESCRIPTION',
     $        'Total '//trim(KSDFT)//' energy')

      if (iUHF.eq.0) then

*     RHF *
***********
        call mh5_init_attr (wfn_fileid,
     $          'ORBITAL_TYPE', trim(KSDFT)//'-RHF')

*     typestring
        wfn_tpidx = mh5_create_dset_str(wfn_fileid,
     $          'MO_TYPEINDICES', 1, [nnB], 1)
        call mh5_init_attr(wfn_tpidx, 'DESCRIPTION',
     $          'Type index of the molecular orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital coefficients
        wfn_mocoef = mh5_create_dset_real(wfn_fileid,
     $          'MO_VECTORS', 1, [nBB])
        call mh5_init_attr(wfn_mocoef, 'DESCRIPTION',
     $          'Coefficients of the molecular orbitals, '//
     $          'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
*     molecular orbital occupation numbers
        wfn_occnum = mh5_create_dset_real(wfn_fileid,
     $          'MO_OCCUPATIONS', 1, [nnB])
        call mh5_init_attr(wfn_occnum, 'DESCRIPTION',
     $          'Occupation numbers of the molecular orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital energies
        wfn_orbene = mh5_create_dset_real(wfn_fileid,
     $          'MO_ENERGIES', 1, [nnB])
        call mh5_init_attr(wfn_orbene, 'DESCRIPTION',
     $          'Orbital energies of the molecular orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')

      else

*     UHF *
***********
        call mh5_init_attr (wfn_fileid,
     $          'ORBITAL_TYPE', trim(KSDFT)//'-UHF')

*     typestring
        wfn_tpidx = mh5_create_dset_str(wfn_fileid,
     $          'MO_TYPEINDICES', 1, [nnB], 1)
        call mh5_init_attr(wfn_tpidx, 'DESCRIPTION',
     $          'Type index of the natural orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital coefficients
        wfn_mocoef = mh5_create_dset_real(wfn_fileid,
     $          'MO_VECTORS', 1, [nBB])
        call mh5_init_attr(wfn_mocoef, 'DESCRIPTION',
     $          'Coefficients of the natural orbitals, '//
     $          'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
*     molecular orbital occupation numbers
        wfn_occnum = mh5_create_dset_real(wfn_fileid,
     $          'MO_OCCUPATIONS', 1, [nnB])
        call mh5_init_attr(wfn_occnum, 'DESCRIPTION',
     $          'Occupation numbers of the natural orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital energies
        wfn_orbene = mh5_create_dset_real(wfn_fileid,
     $          'MO_ENERGIES', 1, [nnB])
        call mh5_init_attr(wfn_orbene, 'DESCRIPTION',
     $          'Orbital energies of the natural orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     typestring: alpha
        wfn_tpidx_a = mh5_create_dset_str(wfn_fileid,
     $          'MO_ALPHA_TYPEINDICES', 1, [nnB], 1)
        call mh5_init_attr(wfn_tpidx_a, 'DESCRIPTION',
     $          'Type index of the alpha orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital coefficients: alpha
        wfn_mocoef_a = mh5_create_dset_real(wfn_fileid,
     $          'MO_ALPHA_VECTORS', 1, [nBB])
        call mh5_init_attr(wfn_mocoef_a, 'DESCRIPTION',
     $          'Coefficients of the alpha orbitals, '//
     $          'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
*     molecular orbital occupation numbers: alpha
        wfn_occnum_a = mh5_create_dset_real(wfn_fileid,
     $          'MO_ALPHA_OCCUPATIONS', 1, [nnB])
        call mh5_init_attr(wfn_occnum_a, 'DESCRIPTION',
     $          'Occupation numbers of the alpha orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital energies: alpha
        wfn_orbene_a = mh5_create_dset_real(wfn_fileid,
     $          'MO_ALPHA_ENERGIES', 1, [nnB])
        call mh5_init_attr(wfn_orbene_a, 'DESCRIPTION',
     $          'Orbital energies of the alpha orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     typestring: beta
        wfn_tpidx_b = mh5_create_dset_str(wfn_fileid,
     $          'MO_BETA_TYPEINDICES', 1, [nnB], 1)
        call mh5_init_attr(wfn_tpidx_b, 'DESCRIPTION',
     $          'Type index of the beta orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital coefficients: beta
        wfn_mocoef_b = mh5_create_dset_real(wfn_fileid,
     $          'MO_BETA_VECTORS', 1, [nBB])
        call mh5_init_attr(wfn_mocoef_b, 'DESCRIPTION',
     $          'Coefficients of the beta orbitals, '//
     $          'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
*     molecular orbital occupation numbers: beta
        wfn_occnum_b = mh5_create_dset_real(wfn_fileid,
     $          'MO_BETA_OCCUPATIONS', 1, [nnB])
        call mh5_init_attr(wfn_occnum_b, 'DESCRIPTION',
     $          'Occupation numbers of the beta orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
*     molecular orbital energies: beta
        wfn_orbene_b = mh5_create_dset_real(wfn_fileid,
     $          'MO_BETA_ENERGIES', 1, [nnB])
        call mh5_init_attr(wfn_orbene_b, 'DESCRIPTION',
     $          'Orbital energies of the beta orbitals '//
     $          'arranged as blocks of size [NBAS(i)], i=1,#irreps')
      end if

#endif
      end
