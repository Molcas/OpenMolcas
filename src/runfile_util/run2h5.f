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
#ifdef _HDF5_
      subroutine run2h5_molinfo(fileid)
*     SVC: read basic molecular information from the RunFile
*     and 1-electron integral file and write it to the HDF5
*     file specified with fileid.
*
*     Attributes:
*       NSYM, IRREP_LABELS, POTNUC, NBAS, NATOMS_UNIQUE, NATOMS_ALL,
*       NPRIM
*     Datasets:
*       CENTER_LABELS, CENTER_CHARGES, CENTER_COORDINATES,
*       BASIS_FUNCTION_IDS, DESYM_CENTER_LABELS, DESYM_CENTER_CHARGES,
*       DESYM_CENTER_COORDINATES, DESYM_BASIS_FUNCTION_IDS,
*       DESYM_MATRIX, PRIMITIVE_IDS, PRIMITIVES

      implicit none
      integer :: fileid
#  include "Molcas.fh"
#  include "stdalloc.fh"
#  include "mh5.fh"

      integer :: natoms
      real*8, allocatable :: charges(:), coord(:,:)
      integer, allocatable :: atnums(:)

      integer :: nsym, isym
      character :: lIrrep(24)

      integer :: nbas(8)
      integer :: nb, nbast, nbast1, nbast2

      real*8 :: potnuc

      character(LENIN), allocatable :: atomlbl(:)
      character(LENIN4), allocatable :: desym_atomlbl(:)

      integer, allocatable :: basis_ids(:,:)
      integer, allocatable :: desym_basis_ids(:,:)

      integer :: mcentr

      real*8, allocatable :: desym_matrix(:)

      integer :: nPrim
      integer, allocatable :: PrimIDs(:,:)
      real*8, allocatable :: Primitives(:,:)

      integer :: i,j
      integer, allocatable :: QMMap(:)

      integer :: dsetid

*     symmetry information
      call get_iScalar('nSym',nSym)
      call mh5_init_attr (fileid,'NSYM', nSym)
      Call Get_cArray('Irreps',lIrrep,24)
      call mh5_init_attr (fileid,'IRREP_LABELS', 1, [nSym], lIrrep, 3)

*     orbital partitions
      Call Get_iArray('nBas',nBas,nSym)
      call mh5_init_attr (fileid,'NBAS', 1, [nSym], nBas)

      nbast=0
      nbast1=0
      nbast2=0
      do isym=1,nsym
        nb=nbas(isym)
        nbast=nbast+nb
        nbast1=nbast1+(nb*(nb+1))/2
        nbast2=nbast2+nb**2
      end do

*     nuclear potential
      Call Get_dScalar('PotNuc',PotNuc)
      call mh5_init_attr (fileid,'POTNUC', PotNuc)

      Call Get_iScalar('Unique centers',nAtoms)
      call mh5_init_attr (fileid,'NATOMS_UNIQUE', nAtoms)

*     atom labels
      dsetid = mh5_create_dset_str(fileid,
     $        'CENTER_LABELS', 1, [nAtoms], LENIN)
      call mh5_init_attr(dsetid, 'description',
     $          'Unique center labels '//
     $          'arranged as one [NATOMS_UNIQUE] block')
      call mma_allocate(atomlbl,nAtoms)
      Call Get_cArray('Un_cen Names',atomlbl,LENIN*nAtoms)
      call mh5_put_dset(dsetid,atomlbl)
      call mma_deallocate(atomlbl)
      call mh5_close_dset(dsetid)

*     atom numbers
      dsetid = mh5_create_dset_int(fileid,
     $        'CENTER_ATNUMS', 1, [nAtoms])
      call mh5_init_attr(dsetid, 'description',
     $        'Atomic numbers, storead as '//
     $        'array of size [NATOMS_UNIQUE]')
      call mma_allocate(atnums,nAtoms)
      Call Get_iArray('Un_cen Charge',atnums,nAtoms)
      call mh5_put_dset(dsetid,atnums)
      call mma_deallocate(atnums)
      call mh5_close_dset(dsetid)

*     atom charges
      dsetid = mh5_create_dset_real(fileid,
     $        'CENTER_CHARGES', 1, [nAtoms])
      call mh5_init_attr(dsetid, 'description',
     $        'Nuclear charges, storead as '//
     $        'array of size [NATOMS_UNIQUE]')
      call mma_allocate(charges,nAtoms)
      Call Get_dArray('Un_cen Effective Charge',charges,nAtoms)
      call mh5_put_dset(dsetid,charges)
      call mma_deallocate(charges)
      call mh5_close_dset(dsetid)

*     atom coordinates
      dsetid = mh5_create_dset_real(fileid,
     $        'CENTER_COORDINATES', 2, [3,nAtoms])
      call mh5_init_attr(dsetid, 'description',
     $        'Atom coordinates, matrix of size [NATOMS_UNIQUE,3], '//
     $        'stored with atom index varying slowest')
      call mma_allocate(coord,3,nAtoms)
      Call Get_dArray('Un_cen Coordinates',coord,3*nAtoms)
      call mh5_put_dset_array_real(dsetid,coord)
      call mma_deallocate(coord)
      call mh5_close_dset(dsetid)

*     unique basis names
      dsetid = mh5_create_dset_int(fileid,
     $        'BASIS_FUNCTION_IDS', 2, [4,nbast])
      call mh5_init_attr(dsetid, 'description',
     $        'Unique basis function IDs (c,n,l,m) '//
     $        'arranged as blocks of size [4*NBAS(i)], i=1,#irreps')
      call mma_allocate(basis_ids,4,nbast)
      Call Get_iArray('Basis IDs',basis_ids,4*nbast)
      call mh5_put_dset_array_int(dsetid,basis_ids)
      call mma_deallocate(basis_ids)
      call mh5_close_dset(dsetid)

      Call get_iScalar('LP_nCenter',mCentr)

      IF (NSYM.GT.1) THEN

        call mh5_init_attr (fileid,'NATOMS_ALL', mCentr)

*     desymmetrized atom labels
        dsetid = mh5_create_dset_str(fileid,
     $          'DESYM_CENTER_LABELS', 1, [mcentr], LENIN4)
        call mh5_init_attr(dsetid, 'description',
     $          'Desymmetrized center labels '//
     $          'arranged as one [NATOMS_ALL] block')
        call mma_allocate(desym_atomlbl,mcentr)
        Call get_cArray('LP_L',desym_atomlbl,(LENIN4)*mCentr)
        call mh5_put_dset(dsetid,desym_atomlbl)
        call mma_deallocate(desym_atomlbl)
        call mh5_close_dset(dsetid)

*     desymmetrized atom numbers
        dsetid = mh5_create_dset_int(fileid,
     $          'DESYM_CENTER_ATNUMS', 1, [MCENTR])
        call mh5_init_attr(dsetid, 'description',
     $          'Desymmetrized atomic numbers, '//
     $          'stored as array of size [NATOMS_ALL]')
        call mma_allocate(atnums,mcentr)
        call get_iArray('LP_A',atnums,mcentr)
        call mh5_put_dset(dsetid,atnums)
        call mma_deallocate(atnums)
        call mh5_close_dset(dsetid)

*     desymmetrized atom charges
        dsetid = mh5_create_dset_real(fileid,
     $          'DESYM_CENTER_CHARGES', 1, [MCENTR])
        call mh5_init_attr(dsetid, 'description',
     $          'Desymmetrized center charges, '//
     $          'stored as array of size [NATOMS_ALL]')
        call mma_allocate(charges,mcentr)
        call get_dArray('LP_Q',charges,mcentr)
        call mh5_put_dset(dsetid,charges)
        call mma_deallocate(charges)
        call mh5_close_dset(dsetid)

*     desymmetrized atom coordinates
        dsetid = mh5_create_dset_real(fileid,
     $          'DESYM_CENTER_COORDINATES', 2, [3,MCENTR])
        call mh5_init_attr(dsetid, 'description',
     $          'Desymmetrized coordinates, size [NATOMS_ALL,3], '//
     $          'stored with atom index varying slowest')
        call mma_allocate(coord,3,MCENTR)
        call get_dArray('LP_Coor',coord,3*mcentr)
        call mh5_put_dset_array_real(dsetid,coord)
        call mma_deallocate(coord)
        call mh5_close_dset(dsetid)

        dsetid = mh5_create_dset_int(fileid,
     $          'DESYM_BASIS_FUNCTION_IDS', 2, [4,nBast])
        call mh5_init_attr(dsetid, 'description',
     $          'Basis function IDs (desymmetrized) (c,n,l,m) '//
     $          'arranged as one [4*NBAST] block, NBAST=sum(NBAS)')
        call mma_allocate(desym_basis_ids, 4, nbast)
        Call get_iArray('Desym Basis IDs',desym_basis_ids,4*nbast)
        call mh5_put_dset_array_int(dsetid,desym_basis_ids)
        call mma_deallocate(desym_basis_ids)
        call mh5_close_dset(dsetid)

*     basis function conversion matrix
        dsetid = mh5_create_dset_real(fileid,
     $          'DESYM_MATRIX', 1, [nbast**2])
        call mh5_init_attr(dsetid, 'description',
     $          'Symmetrization matrix for the basis functions '//
     $          'arranged as a [NBAST,NBAST] block, NBAST=sum(NBAS), '//
     $          'fast index corresponds to desymmetrized basis.')
        call mma_allocate(desym_matrix, nbast**2)
        Call get_dArray('SM',desym_matrix,nbast**2)
        call mh5_put_dset(dsetid,desym_matrix)
        call mma_deallocate(desym_matrix)
        call mh5_close_dset(dsetid)

      END IF

      call get_iScalar('nPrim',nPrim)
      call mh5_init_attr (fileid,'NPRIM', NPRIM)

*     radial primitive identifications
      dsetid = mh5_create_dset_int(fileid,
     $        'PRIMITIVE_IDS', 2, [3,nPrim])
      call mh5_init_attr(dsetid, 'description',
     $        'Primitive IDs, arranged as an '//
     $        'array of size [3*NPRIM], with consecutive '//
     $        'center_id, angmom, shell_id (C1 2s <-> 1,0,2)')
      call mma_allocate(PrimIDs,3,nPrim)
      call get_iArray('primitive ids', PrimIDs, 3*nPrim)
*     only QM atoms have basis functions,
*     remap the center_id numbers to skip MM atoms
      call mma_allocate(QMMap,mcentr)
      call get_iArray('IsMM Atoms', QMMap, mcentr)
      j=0
      do i=1,mcentr
        if (QMMap(i).eq.0) then
          j = j+1
          QMMap(j) = i
        end if
      end do
      if (j.lt.mcentr) then
        do i=1,nPrim
          PrimIDs(1,i) = QMMap(PrimIDs(1,i))
        end do
      end if
      call mma_deallocate(QMMap)
      call mh5_put_dset_array_int(dsetid,PrimIDs)
      call mma_deallocate(PrimIDs)
      call mh5_close_dset(dsetid)

*     atom maximum angular momenta
      dsetid = mh5_create_dset_real(fileid,
     $        'PRIMITIVES', 2, [2,nPrim])
      call mh5_init_attr(dsetid, 'description',
     $        'Primitives, arranged as an '//
     $        'array of size [2*NPRIM], with consecutive '//
     $        'exponent, contraction coefficient')
      call mma_allocate(Primitives,2,nPrim)
      call get_dArray('primitives', Primitives, 2*nPrim)
      call mh5_put_dset_array_real(dsetid,Primitives)
      call mma_deallocate(Primitives)
      call mh5_close_dset(dsetid)

      end
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_run2h5_molinfo()
      End
#endif
