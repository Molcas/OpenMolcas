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
      subroutine cre_dyn
*     MV: Create a dynamix file. If another .dyn.h5 file already
*     exists, it will be overwritten.
      implicit none
#ifdef _HDF5_
#  include "Molcas.fh"
#  include "dyn.fh"
#  include "stdalloc.fh"
      integer :: natoms,nh,nsym,nstates,nconfs,
     $           dyn_dsetid,surf_dsetid,wfn_fileid,ii
      real*8, allocatable :: coord(:,:), ener(:), ciarray(:)
      character(len=LENIN), allocatable :: atomlbl(:)
      logical :: found
      complex*16, allocatable :: amatrix(:)

*     Create a new dynamix file!
      dyn_fileid = mh5_create_file('DYN')

*     Set module type
      call mh5_init_attr (dyn_fileid,'MOLCAS_MODULE', 'DYNAMIX')

*     Copy basic molecular information to the HDF5 file
*     Morgane Vacher: Need to be able to modify the geom
*     but could not make it work with run2h5
*      call run2h5_molinfo(dyn_fileid)
*     Symmetry
      call mh5_init_attr (dyn_fileid,'NSYM', 1)
      call get_iScalar('nSym',nsym)
      if (nsym .gt. 1) then
        Call get_iScalar('LP_nCenter',natoms)
      else
        Call Get_iScalar('Unique atoms',natoms)
      endif
*     Number of atoms
      call mh5_init_attr (dyn_fileid,'NATOMS_UNIQUE', natoms)
*     Nuclear geometry
      dyn_geom = mh5_create_dset_real(dyn_fileid,
     $          'CENTER_COORDINATES', 1, [3*natoms])
      call mh5_init_attr(dyn_geom, 'DESCRIPTION',
     $          'Geometry in cartesians coordinates '//
     $          'at the current time step')
      call mma_allocate(coord,3,natoms)
      if (nsym .gt. 1) then
        call get_dArray('LP_Coor',coord,3*natoms)
      else
        Call Get_dArray('Unique Coordinates',coord,3*natoms)
      endif
      call mh5_put_dset_array_real(dyn_geom,coord)
      call mma_deallocate(coord)
*     Atom labels
      dyn_dsetid = mh5_create_dset_str(dyn_fileid,
     $          'CENTER_LABELS', 1, [natoms], LENIN)
      call mh5_init_attr(dyn_dsetid, 'DESCRIPTION',
     $          'Center labels '//
     $          'arranged as a [NATOMS] block')
      call mma_allocate(atomlbl,natoms)
      if (nsym .gt. 1) then
        call get_cArray('LP_L',atomlbl,LENIN4*natoms)
      else
        Call Get_cArray('Unique Atom Names',atomlbl,LENIN*natoms)
      endif
      call mh5_put_dset(dyn_dsetid,atomlbl)
      call mma_deallocate(atomlbl)
      call mh5_close_dset(dyn_dsetid)

*     Current time
      dyn_time = mh5_create_dset_real (dyn_fileid,'TIME')
      call mh5_init_attr(dyn_time, 'DESCRIPTION',
     $        'Current time of the Molecular Dynamics')
*     Time step
      dyn_dt = mh5_create_dset_real (dyn_fileid,'TIME_STEP')
      call mh5_init_attr(dyn_dt, 'DESCRIPTION',
     $        'Time step of the Molecular Dynamics')

*     Total energy at t=0
      dyn_etot0 = mh5_create_dset_real (dyn_fileid,'ETOT_0')
      call mh5_init_attr(dyn_etot0, 'DESCRIPTION',
     $        'Total energy at t=0')
*     Total energy at current time
      dyn_etot = mh5_create_dset_real (dyn_fileid,'ETOT')
      call mh5_init_attr(dyn_etot, 'DESCRIPTION',
     $        'Total energy at current time')

*     Velocities
      dyn_vel = mh5_create_dset_real(dyn_fileid,
     $          'VELOCITIES', 1, [3*natoms])
      call mh5_init_attr(dyn_vel, 'DESCRIPTION',
     $          'Velocities in cartesians coordinates '//
     $          'at the current time step')

*     NoseHoover
      nh=6
      dyn_nh = mh5_create_dset_real(dyn_fileid,
     $          'NOSEHOOVER', 1, [nh])
      call mh5_init_attr(dyn_nh, 'DESCRIPTION',
     $          'NoseHoover degrees of freedom')

*     MaxHop
*     Morgane Vacher: Dataset only created if needed since its existence serves as a flag.
      call qpg_iScalar('MaxHops',Found)
      if (Found) then
         call get_iScalar('MaxHops',ii)
         call mh5_init_dset(dyn_fileid,'MAX_HOP',ii)
      end if

*     Isotopes
      dyn_mass = mh5_create_dset_real(dyn_fileid,
     $          'MASSES', 1, [natoms])
      call mh5_init_attr(dyn_mass, 'DESCRIPTION',
     $          'Atomic masses, in a.u.')

*     The following variables are relevant to the SURFACEHOP module.
*     They are not modified by the DYNAMIX module but saved here on the H5 file.
*     This way, the informations saved are consistent, ie belong to the same step.

*     Seed number
      call qpg_iscalar('Seed',Found)
      if (Found) then
        surf_dsetid = mh5_create_dset_int (dyn_fileid,'SEED')
        call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'Seed number')
        call get_iscalar('Seed',ii)
        call mh5_put_dset(surf_dsetid,ii)
        call mh5_close_dset(surf_dsetid)
      endif

*     Number of hops
      call qpg_iscalar('Number of Hops',Found)
      if (Found) then
        surf_dsetid = mh5_create_dset_int (dyn_fileid,'NO. OF HOPS')
        call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'Number of hops')
        call get_iscalar('Number of Hops',ii)
        call mh5_put_dset(surf_dsetid,ii)
        call mh5_close_dset(surf_dsetid)
      endif

*     Maximum number of hops in Tully
      call qpg_iscalar('MaxHopsTully',Found)
      if (Found) then
        surf_dsetid = mh5_create_dset_int (dyn_fileid,'MAX_HOP_TULLY')
        call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'Maximum number of hops in Tully algorithm')
        call get_iscalar('MaxHopsTully',ii)
        call mh5_put_dset(surf_dsetid,ii)
        call mh5_close_dset(surf_dsetid)
      endif

*     Relax CASSCF root
      call qpg_iscalar('Relax CASSCF root',Found)
      if (Found) then
        surf_dsetid = mh5_create_dset_int (dyn_fileid,'RELAX CAS ROOT')
        call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'Relax CASSCF root')
        call get_iscalar('Relax CASSCF root',ii)
        call mh5_put_dset(surf_dsetid,ii)
        call mh5_close_dset(surf_dsetid)
      endif

*     Read number of states and configurations from rasscf.h5 file
      wfn_fileid = mh5_open_file_r('RASWFN')
      if (mh5_exists_attr(wfn_fileid,'NSTATES')
     $  .and. mh5_exists_attr(wfn_fileid,'NCONF')) then
        call mh5_fetch_attr(wfn_fileid,'NSTATES',nstates)
        call mh5_fetch_attr(wfn_fileid,'NCONF',nconfs)
        call mh5_init_attr(dyn_fileid,'NSTATES',nstates)
        call mh5_init_attr(dyn_fileid,'NCONFS',nconfs)

*     Energies at the previous step
        call qpg_darray('VenergyP',Found,nstates)
        if (Found) then
          surf_dsetid = mh5_create_dset_real (dyn_fileid,'ENERG PREV',
     $            1, [nstates])
          call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'Potential energies at the previous time step')
          call mma_allocate(ener,nstates)
          call get_darray('VenergyP',ener,nstates)
          call mh5_put_dset(surf_dsetid,ener)
          call mma_deallocate(ener)
          call mh5_close_dset(surf_dsetid)
        endif

*     CI coeffs at the previous step
        call qpg_darray('AllCIP',Found,nstates*nconfs)
        if (Found) then
          surf_dsetid = mh5_create_dset_real (dyn_fileid,'CI PREV',
     $            1, [nstates*nconfs])
          call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'CI coeffs at the previous time step')
          call mma_allocate(ciarray,nstates*nconfs)
          call get_darray('AllCIP',ciarray,nstates*nconfs)
          call mh5_put_dset(surf_dsetid,ciarray)
          call mma_deallocate(ciarray)
          call mh5_close_dset(surf_dsetid)
        endif

*     CI coeffs at the step before the previous step
        call qpg_darray('AllCIPP',Found,nstates*nconfs)
        if (Found) then
          surf_dsetid = mh5_create_dset_real (dyn_fileid,'CI PPREV',
     $            1, [nstates*nconfs])
          call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'CI coeffs at the step before the previous time step')
          call mma_allocate(ciarray,nstates*nconfs)
          call get_darray('AllCIPP',ciarray,nstates*nconfs)
          call mh5_put_dset(surf_dsetid,ciarray)
          call mma_deallocate(ciarray)
          call mh5_close_dset(surf_dsetid)
        endif

*     A matrix V
        call Qpg_zArray('AmatrixV', Found, nstates*nstates)
        if (Found) then
          call mma_allocate(Amatrix,nstates*nstates)
          call get_zarray('AmatrixV',Amatrix,NSTATES*NSTATES)
*     HDF5 format does not deal with complex numbers so split manually into
*     real and imaginary parts and save 2 datasets
*         Real part
          surf_dsetid = mh5_create_dset_real (dyn_fileid,'AMATRIXV-R',
     $            1, [nstates*nstates])
          call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'real part of AmatrixV')
          call mh5_put_dset(surf_dsetid,REAL(Amatrix))
          call mh5_close_dset(surf_dsetid)
*         Imaginary part
          surf_dsetid = mh5_create_dset_real (dyn_fileid,'AMATRIXV-I',
     $            1, [nstates*nstates])
          call mh5_init_attr(surf_dsetid, 'DESCRIPTION',
     $        'imaginary part of AmatrixV')
          call mh5_put_dset(surf_dsetid,AIMAG(Amatrix))
          call mh5_close_dset(surf_dsetid)
          call mma_deallocate(Amatrix)
        endif


      endif
      call mh5_close_file(wfn_fileid)

#endif
      end
