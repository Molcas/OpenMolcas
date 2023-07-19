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
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

module Chkpnt

#ifdef _HDF5_
use mh5, only: mh5_close_attr, mh5_close_dset, mh5_close_file, mh5_create_attr_int, mh5_create_dset_int, mh5_create_dset_real, &
               mh5_create_dset_str, mh5_create_file, mh5_fetch_attr, mh5_get_attr, mh5_init_attr, mh5_is_hdf5, mh5_open_attr, &
               mh5_open_dset, mh5_open_file_rw, mh5_put_attr, mh5_put_dset, mh5_resize_dset
use Definitions, only: wp
#endif
use Definitions, only: iwp

implicit none
private

# ifdef _HDF5_
character(len=9), parameter :: basename = 'SLAPAFCHK'

integer(kind=iwp) :: chkpnt_coor, chkpnt_ener, chkpnt_force, chkpnt_hess, chkpnt_id, chkpnt_iter, chkpnt_new, Iter_all
character(len=12) :: filename
#endif

public :: Chkpnt_close, Chkpnt_open, Chkpnt_update, Chkpnt_update_MEP

contains

subroutine Chkpnt_open()
# ifdef _HDF5_
  use Symmetry_Info, only: nIrrep
  use Slapaf_Info, only: Coor, IRC, iter
  integer(kind=iwp) :: tmp
  logical(kind=iwp) :: create
  character(len=3) :: level

  Iter_all = Iter
  create = .true.
  call GetEnvF('EMIL_InLoop',level)
  if ((level == '0') .or. (level == '1')) level = ''
  filename = basename//trim(level)

  if (((Iter > 1) .or. (IRC == -1)) .and. mh5_is_hdf5(filename)) then
    create = .false.
    chkpnt_id = mh5_open_file_rw(filename)
    chkpnt_iter = mh5_open_attr(chkpnt_id,'ITERATIONS')
    chkpnt_ener = mh5_open_dset(chkpnt_id,'ENERGIES')
    chkpnt_coor = mh5_open_dset(chkpnt_id,'COORDINATES')
    chkpnt_new = mh5_open_dset(chkpnt_id,'CENTER_COORDINATES')
    chkpnt_force = mh5_open_dset(chkpnt_id,'FORCES')
    chkpnt_hess = mh5_open_dset(chkpnt_id,'HESSIAN')
    call mh5_fetch_attr(chkpnt_id,'NSYM',tmp)
    if (tmp /= nIrrep) create = .true.
    call mh5_fetch_attr(chkpnt_id,'NATOMS_UNIQUE',tmp)
    if (tmp /= size(Coor,2)) create = .true.
    call mh5_fetch_attr(chkpnt_id,'ITERATIONS',tmp)
    if (IRC == -1) then
      if (Iter >= tmp) create = .true.
      Iter_all = tmp+1
    else
      if (tmp >= Iter) create = .true.
    end if
  end if
  if (create) then
    if ((Iter > 1) .or. (IRC == -1)) then
      call WarningMessage(2,'The HDF5 file does not exist or is inconsistent')
      call AbEnd()
    else
      call Chkpnt_init()
    end if
  end if
# endif
end subroutine Chkpnt_open

#ifdef _HDF5_
subroutine Chkpnt_init()
  use Phase_Info, only: iPhase
  use Symmetry_Info, only: nIrrep
  use Index_Functions, only: nTri_Elem
  use Slapaf_Info, only: AtomLbl, Coor, dMass, dMEPStep, iCoSet, MEP, nDimBC, nStab, rMEP, Smmtrc
  use stdalloc, only: mma_allocate, mma_deallocate
# include "Molcas.fh"
  character :: lIrrep(24)
  integer(kind=iwp) :: dsetid, i, j, k, mAtom
  integer(kind=iwp), allocatable :: desym(:,:), symdof(:,:)
  real(kind=wp), allocatable :: charges(:)

  chkpnt_id = mh5_create_file(filename)

  call mh5_init_attr(chkpnt_id,'MOLCAS_MODULE','SLAPAF')

  ! symmetry information
  call mh5_init_attr(chkpnt_id,'NSYM',nIrrep)
  call Get_cArray('Irreps',lIrrep,24)
  call mh5_init_attr(chkpnt_id,'IRREP_LABELS',1,[nIrrep],lIrrep,3)

  call mh5_init_attr(chkpnt_id,'NATOMS_UNIQUE',size(Coor,2))

  call mh5_init_attr(chkpnt_id,'DOF',nDimBC)

  ! atom labels
  dsetid = mh5_create_dset_str(chkpnt_id,'CENTER_LABELS',1,[size(Coor,2)],LenIn)
  call mh5_init_attr(dsetid,'DESCRIPTION','Unique center labels arranged as one [NATOMS_UNIQUE] block')
  call mh5_put_dset(dsetid,AtomLbl)
  call mh5_close_dset(dsetid)

  ! atom masses
  dsetid = mh5_create_dset_real(chkpnt_id,'CENTER_MASSES',1,[size(Coor,2)])
  call mh5_init_attr(dsetid,'DESCRIPTION','Nuclear masses, stored as array of size [NATOMS_UNIQUE]')
  call mh5_put_dset(dsetid,dMass)
  call mh5_close_dset(dsetid)

  ! atom charges
  dsetid = mh5_create_dset_real(chkpnt_id,'CENTER_CHARGES',1,[size(Coor,2)])
  call mh5_init_attr(dsetid,'DESCRIPTION','Nuclear charges, stored as array of size [NATOMS_UNIQUE]')
  call mma_allocate(charges,size(Coor,2))
  call Get_dArray('Nuclear Charge',charges,size(Coor,2))
  call mh5_put_dset(dsetid,charges)
  call mma_deallocate(charges)
  call mh5_close_dset(dsetid)

  ! number of iterations
  chkpnt_iter = mh5_create_attr_int(chkpnt_id,'ITERATIONS')

  ! atom coordinates (new iteration)
  !   use the same dataset name as in run2hdf5
  chkpnt_new = mh5_create_dset_real(chkpnt_id,'CENTER_COORDINATES',2,[3,size(Coor,2)])
  call mh5_init_attr(chkpnt_new,'DESCRIPTION','Atom coordinates for new iteration, matrix of size [NATOMS_UNIQUE,3], stored '// &
                     'with atom index varying slowest')

  if (nIrrep > 1) then

    mAtom = 0
    do i=1,size(Coor,2)
      mAtom = mAtom+nIrrep/nStab(i)
    end do
    call mma_allocate(desym,4,mAtom)
    call mma_allocate(symdof,2,nDimBC)
    mAtom = 0
    k = 0
    do i=1,size(Coor,2)
      do j=0,nIrrep/nStab(i)-1
        mAtom = mAtom+1
        desym(1,mAtom) = i
        desym(2,mAtom) = iPhase(1,iCoSet(j,i))
        desym(3,mAtom) = iPhase(2,iCoSet(j,i))
        desym(4,mAtom) = iPhase(3,iCoSet(j,i))
      end do
      do j=1,3
        if (.not. Smmtrc(j,i)) cycle
        k = k+1
        symdof(1,k) = i
        symdof(2,k) = j
      end do
    end do

    ! total number of atoms
    call mh5_init_attr(chkpnt_id,'NATOMS_ALL',mAtom)

    ! desymmetrization factors
    dsetid = mh5_create_dset_int(chkpnt_id,'DESYM_FACTORS',2,[4,mAtom])
    call mh5_init_attr(dsetid,'DESCRIPTION','Factors for obtaining all coordinates, matrix of size [NATOMS_ALL,4], each row '// &
                       'contains the unique atom index and the factors with which to multiply the x,y,z coordinates')
    call mh5_put_dset(dsetid,desym)
    call mh5_close_dset(dsetid)
    call mma_deallocate(desym)

    ! symmetry-unique degrees of freedom (Cartesian indices)
    dsetid = mh5_create_dset_int(chkpnt_id,'DOF_INDICES',2,[2,nDimBC])
    call mh5_init_attr(dsetid,'DESCRIPTION','Indices of the Cartesian degrees of freedom, matrix of size [DOF, 2], each row '// &
                       'contains the atom index and the Cartesian index (1=x, 2=y, 3=z)')
    call mh5_put_dset(dsetid,symdof)
    call mh5_close_dset(dsetid)
    call mma_deallocate(symdof)
  end if

  ! iteration data:

  ! energies
  chkpnt_ener = mh5_create_dset_real(chkpnt_id,'ENERGIES',1,[0],dyn=.true.)
  call mh5_init_attr(chkpnt_ener,'DESCRIPTION','Energies for all iterations as a matrix of size [ITERATIONS]')

  ! atom coordinates
  chkpnt_coor = mh5_create_dset_real(chkpnt_id,'COORDINATES',3,[3,size(Coor,2),0],dyn=.true.)
  call mh5_init_attr(chkpnt_coor,'DESCRIPTION','Atom coordinates, matrix of size [ITERATIONS,NATOMS_UNIQUE,3], stored with '// &
                     'iteration varying slowest, then atom index')

  ! Cartesian forces (F = -g)
  chkpnt_force = mh5_create_dset_real(chkpnt_id,'FORCES',3,[3,size(Coor,2),0],dyn=.true.)
  call mh5_init_attr(chkpnt_force,'DESCRIPTION','Cartesian forces, matrix of size [ITERATIONS,NATOMS_UNIQUE,3], stored with '// &
                     'iteration varying slowest, then atom index')

  ! Cartesian Hessian
  chkpnt_hess = mh5_create_dset_real(chkpnt_id,'HESSIAN',1,[nTri_Elem(nDimBC)])
  call mh5_init_attr(chkpnt_hess,'DESCRIPTION','Cartesian Hessian in triangular form, as a vector of size [DOF*(DOF+1)/2]')

  ! MEP/IRC information
  if (MEP .or. rMEP) then
    call mh5_init_attr(chkpnt_id,'MEP_STEP',dMEPStep)
    call mh5_init_attr(chkpnt_id,'MEP_ITERATIONS',0)

    dsetid = mh5_create_dset_int(chkpnt_id,'MEP_INDICES',1,[0],dyn=.true.)
    call mh5_init_attr(dsetid,'DESCRIPTION','Iteration number for each converged MEP step, as a vector of size [MEP_ITERATIONS]')
  end if
end subroutine Chkpnt_init
#endif

subroutine Chkpnt_update()
# ifdef _HDF5_
  use Slapaf_Info, only: Cx, Energy, Gx, iter, nDimBC
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp) :: i, ij, j
  logical(kind=iwp) :: Found
  real(kind=wp), allocatable :: Hss_X(:)

  call Qpg_dArray('Hss_X',Found,i)
  if (Found) then
    if (i /= nDimBC**2) then
      call WarningMessage(2,'Hessian with wrong dimension')
      call AbEnd()
    end if
    call mma_allocate(Hss_X,i)
    call Get_dArray('Hss_X',Hss_X,i)
    ij = 0
    do i=1,nDimBC
      do j=1,i
        ij = ij+1
        Hss_X(ij) = Hss_X(nDimBC*(i-1)+j)
      end do
    end do
  end if

  ! iterations
  call mh5_put_attr(chkpnt_iter,Iter_all)
  ! energies
  call mh5_resize_dset(chkpnt_ener,[Iter_all])
  call mh5_put_dset(chkpnt_ener,Energy(Iter:Iter),[1],[Iter_all-1])
  ! coordinates
  call mh5_resize_dset(chkpnt_coor,[3,size(Cx,2),Iter_all])
  call mh5_put_dset(chkpnt_coor,Cx(:,:,Iter),[3,size(Cx,2),1],[0,0,Iter_all-1])
  ! new coordinates
  call mh5_put_dset(chkpnt_new,Cx(1,1,Iter+1))
  ! forces
  call mh5_resize_dset(chkpnt_force,[3,size(Cx,2),Iter_all])
  call mh5_put_dset(chkpnt_force,Gx(:,:,Iter),[3,size(Cx,2),1],[0,0,Iter_all-1])
  ! Hessian
  if (Found) then
    call mh5_put_dset(chkpnt_hess,Hss_X(1))
    call mma_deallocate(Hss_X)
  end if
# endif
end subroutine Chkpnt_update

subroutine Chkpnt_update_MEP(SaveMEP,IRCRestart)
  logical(kind=iwp), intent(in) :: SaveMEP, IRCRestart
# ifdef _HDF5_
  integer(kind=iwp) :: attrid, dsetid, iMEP

  if (IRCRestart) call mh5_init_attr(chkpnt_id,'IRC_RESTART',Iter_all+1)
  if (SaveMEP) then
    attrid = mh5_open_attr(chkpnt_id,'MEP_ITERATIONS')
    call mh5_get_attr(attrid,iMEP)
    iMEP = iMEP+1
    call mh5_put_attr(attrid,iMEP)
    call mh5_close_attr(attrid)
    dsetid = mh5_open_dset(chkpnt_id,'MEP_INDICES')
    call mh5_resize_dset(dsetid,[iMEP])
    call mh5_put_dset(dsetid,[Iter_all],[1],[iMEP-1])
    call mh5_close_dset(dsetid)
  end if
# else
# include "macros.fh"
  unused_var(SaveMEP)
  unused_var(IRCRestart)
# endif
end subroutine Chkpnt_update_MEP

subroutine Chkpnt_close()
# ifdef _HDF5_
  call mh5_close_file(chkpnt_id)
# endif
end subroutine Chkpnt_close

end module Chkpnt
