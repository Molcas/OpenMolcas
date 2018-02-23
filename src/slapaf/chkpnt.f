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
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
      Module Chkpnt
      Implicit None
#ifdef _HDF5_
#  include "mh5.fh"
#endif
      Character(Len=9) :: basename = 'SLAPAFCHK'
      Character(Len=12) :: filename
      Integer :: chkpnt_id, chkpnt_iter
      Integer :: chkpnt_ener, chkpnt_coor, chkpnt_grad, chkpnt_new
      Integer :: Iter_all

      Contains

      Subroutine Chkpnt_open()
#ifdef _HDF5_
#  include "info_slapaf.fh"
      Character(Len=3) :: level
      Logical :: create
      Integer :: tmp

      Iter_all = Iter
      create = .True.
      Call GetEnvF('EMIL_InLoop',level)
      If ((level.eq.'0').or.(level.eq.'1')) level=''
      filename = basename//Trim(level)

      If (((Iter.gt.1).or.(IRC.eq.-1)).and.mh5_is_hdf5(filename)) Then
        create = .False.
        chkpnt_id = mh5_open_file_rw(filename)
        chkpnt_iter = mh5_open_attr(chkpnt_id, 'ITERATIONS')
        chkpnt_ener = mh5_open_dset(chkpnt_id, 'ENERGIES')
        chkpnt_coor = mh5_open_dset(chkpnt_id, 'COORDINATES')
        chkpnt_new  = mh5_open_dset(chkpnt_id, 'CENTER_COORDINATES')
        chkpnt_grad = mh5_open_dset(chkpnt_id, 'FORCES')
        Call mh5_fetch_attr(chkpnt_id, 'NSYM', tmp)
        If (tmp.ne.nSym) create = .True.
        Call mh5_fetch_attr(chkpnt_id, 'NATOMS_UNIQUE', tmp)
        If (tmp.ne.nsAtom) create = .True.
        Call mh5_fetch_attr(chkpnt_id, 'ITERATIONS', tmp)
        If (IRC.eq.-1) Then
          If (Iter.ge.tmp) create = .True.
          Iter_all = tmp+1
        Else
          If (tmp.ge.Iter) create = .True.
        End If
      End If
      If (create) Then
        If ((Iter.gt.1).or.(IRC.eq.-1)) Then
          Call WarningMessage(2,
     &         'The HDF5 file does not exist or is inconsistent')
          Call AbEnd()
        Else
          Call Chkpnt_init()
        End If
      End If
#endif
      End Subroutine Chkpnt_open

      Subroutine Chkpnt_init()
#ifdef _HDF5_
#  include "info_slapaf.fh"
#  include "WrkSpc.fh"
#  include "stdalloc.fh"
      Character :: lIrrep(24)
      Integer :: dsetid, attrid, mAtom, i, j
      Real*8, Allocatable :: charges(:)
      Integer :: iPhase(3,0:7)
      Data iPhase/ 1, 1, 1,   -1, 1, 1,   1,-1, 1,  -1,-1, 1,
     &             1, 1,-1,   -1, 1,-1,   1,-1,-1,  -1,-1,-1/
      Integer, Allocatable :: desym(:,:)

      chkpnt_id = mh5_create_file(filename)

      Call mh5_init_attr(chkpnt_id, 'MOLCAS_MODULE', 'SLAPAF')

*     symmetry information
      Call mh5_init_attr(chkpnt_id, 'NSYM', nSym)
      Call Get_cArray('Irreps', lIrrep, 24)
      Call mh5_init_attr(chkpnt_id, 'IRREP_LABELS',
     &                   1, [nSym], lIrrep, 3)

      Call mh5_init_attr(chkpnt_id, 'NATOMS_UNIQUE', nsAtom)

*     atom labels
      dsetid = mh5_create_dset_str(chkpnt_id,
     &         'CENTER_LABELS', 1, [nsAtom], LENIN)
      Call mh5_init_attr(dsetid, 'description',
     &     'Unique center labels arranged as one [NATOMS_UNIQUE] block')
      Call mh5_put_dset(dsetid, AtomLbl)
      Call mh5_close_dset(dsetid)

*     atom masses
      dsetid = mh5_create_dset_real(chkpnt_id,
     &         'CENTER_MASSES', 1, [nsAtom])
      Call mh5_init_attr(dsetid, 'description',
     &    'Nuclear masses, stored as array of size [NATOMS_UNIQUE]')
      Call mh5_put_dset(dsetid, Work(ipCM))
      Call mh5_close_dset(dsetid)

*     atom charges
      dsetid = mh5_create_dset_real(chkpnt_id,
     &         'CENTER_CHARGES', 1, [nsAtom])
      Call mh5_init_attr(dsetid, 'description',
     &    'Nuclear charges, stored as array of size [NATOMS_UNIQUE]')
      Call mma_allocate(charges, nsAtom)
      Call Get_dArray('Nuclear Charge', charges, nsAtom)
      Call mh5_put_dset(dsetid, charges)
      Call mma_deallocate(charges)
      Call mh5_close_dset(dsetid)

*     number of iterations
      chkpnt_iter = mh5_create_attr_int(chkpnt_id,
     &              'ITERATIONS')

*     atom coordinates (new iteration)
*       use the same dataset name as in run2hdf5
      chkpnt_new = mh5_create_dset_real(chkpnt_id,
     &             'CENTER_COORDINATES', 2, [3,nsAtom])
      Call mh5_init_attr(chkpnt_new, 'description',
     &     'Atom coordinates for new iteration, matrix of size '//
     &     '[NATOMS_UNIQUE,3], stored with atom index varying slowest')

      If (nSym.gt.1) Then

        mAtom = 0
        Do i=1,nsAtom
          mAtom = mAtom+nSym/nStab(i)
        End Do
        Call mma_allocate(desym, 4, mAtom)
        mAtom = 0
        Do i=1,nsAtom
          Do j=0,nSym/nStab(i)-1
            mAtom = mAtom+1
            desym(1,mAtom) = i
            desym(2,mAtom) = iPhase(1,iCoSet(j,i))
            desym(3,mAtom) = iPhase(2,iCoSet(j,i))
            desym(4,mAtom) = iPhase(3,iCoSet(j,i))
          End Do
        End Do

*     total number of atoms
        Call mh5_init_attr(chkpnt_id, 'NATOMS_ALL', mAtom)

*     desymmetrization factors
        dsetid = mh5_create_dset_int(chkpnt_id,
     &           'DESYM_FACTORS', 2, [4,mAtom])
        Call mh5_init_attr(dsetid, 'description',
     &      'Factors for obtaining all coordinates, matrix of size '//
     &      '[NATOMS_ALL,4], each row contains the unique atom index'//
     &      'and the factors with which to multiply the x,y,z '//
     &      'coordinates')
        Call mh5_put_dset(dsetid, desym(1,1))
        Call mh5_close_dset(dsetid)
        Call mma_deallocate(desym)

      End If

*     iteration data:

*     energies
      chkpnt_ener = mh5_create_dset_real(chkpnt_id,
     &              'ENERGIES', 1, [0], dyn=.True.)
      Call mh5_init_attr(chkpnt_ener, 'description',
     &     'Energies for all iterations as a matrix of size '//
     &     '[ITERATIONS]')

*     atom coordinates
      chkpnt_coor = mh5_create_dset_real(chkpnt_id,
     &              'COORDINATES', 3, [3,nsAtom,0], dyn=.True.)
      Call mh5_init_attr(chkpnt_coor, 'description',
     &     'Atom coordinates, matrix of size [ITERATIONS,'//
     &     'NATOMS_UNIQUE,3], stored with iteration varying slowest, '//
     &     'then atom index')

*     Cartesian forces (F = -g)
      chkpnt_grad = mh5_create_dset_real(chkpnt_id,
     &              'FORCES', 3, [3,nsAtom,0], dyn=.True.)
      Call mh5_init_attr(chkpnt_grad, 'description',
     &     'Cartesian forces, matrix of size [ITERATIONS,'//
     &     'NATOMS_UNIQUE,3], stored with iteration varying slowest, '//
     &     'then atom index')

*     MEP/IRC information
      If (MEP) Then
        Call mh5_init_attr(chkpnt_id, 'MEP_STEP', dMEPStep)
        attrid = mh5_create_attr_int(chkpnt_id, 'MEP_ITERATIONS')

        dsetid = mh5_create_dset_int(chkpnt_id,
     &           'MEP_INDICES', 1, [0], dyn=.True.)
        Call mh5_init_attr(dsetid, 'description',
     &       'Iteration number for each converged MEP step, as a '//
     &       'of size [MEP_ITERATIONS]')
      End If
#endif

      End Subroutine Chkpnt_init

      Subroutine Chkpnt_update()
#ifdef _HDF5_
#  include "info_slapaf.fh"
#  include "WrkSpc.fh"
      Integer :: N3

      N3 = 3*nsAtom
      Call mh5_put_attr(chkpnt_iter, Iter_all)
      Call mh5_resize_dset(chkpnt_ener, [Iter_all])
      Call mh5_put_dset_array_real(chkpnt_ener,
     &     Work(ipEner+(Iter-1)), [1], [Iter_all-1])
      Call mh5_resize_dset(chkpnt_coor, [3,nsAtom,Iter_all])
      Call mh5_put_dset_array_real(chkpnt_coor,
     &     Work(ipCx+N3*(Iter-1)), [3,nsAtom,1], [0,0,Iter_all-1])
      Call mh5_put_dset_array_real(chkpnt_new,
     &     Work(ipCx+N3*Iter), [3,nsAtom], [0,0])
      Call mh5_resize_dset(chkpnt_grad, [3,nsAtom,Iter_all])
      Call mh5_put_dset_array_real(chkpnt_grad,
     &     Work(ipGx+N3*(Iter-1)), [3,nsAtom,1], [0,0,Iter_all-1])
#endif
      End Subroutine Chkpnt_update

      Subroutine Chkpnt_update_MEP(IRCRestart)
      Logical, Intent(In) :: IRCRestart
#ifdef _HDF5_
#  include "info_slapaf.fh"
      Integer :: attrid, dsetid, iMEP

      If (IRCRestart) Then
        Call mh5_init_attr(chkpnt_id, 'IRC_RESTART', Iter_all+1)
      End If
      attrid = mh5_open_attr(chkpnt_id, 'MEP_ITERATIONS')
      Call mh5_get_attr(attrid, iMEP)
      iMEP = iMEP+1
      Call mh5_put_attr(attrid, iMEP)
      Call mh5_close_attr(attrid)
      dsetid = mh5_open_dset(chkpnt_id, 'MEP_INDICES')
      Call mh5_resize_dset(dsetid, [iMEP])
      Call mh5_put_dset_array_int(dsetid, [Iter_all], [1], [iMEP-1])
      Call mh5_close_dset(dsetid)
#endif
      End Subroutine Chkpnt_update_MEP

      Subroutine Chkpnt_close()
#ifdef _HDF5_
      Call mh5_close_file(chkpnt_id)
#endif
      End Subroutine Chkpnt_close

      End Module Chkpnt
