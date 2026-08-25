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
               mh5_create_dset_str, mh5_create_file, mh5_exists_dset, mh5_fetch_attr, mh5_get_attr, mh5_init_attr, mh5_is_hdf5, &
               mh5_open_attr, mh5_open_dset, mh5_open_file_rw, mh5_put_attr, mh5_put_dset, mh5_resize_dset
use Definitions, only: wp
#endif
use Molcas, only: LenIn
use Definitions, only: iwp

implicit none
private

# ifdef _HDF5_
character(len=*), parameter :: basename = 'SLAPAFCHK'

integer(kind=iwp) :: chkpnt_appnadc, chkpnt_coor, chkpnt_ener, chkpnt_force, chkpnt_gd, chkpnt_hess, chkpnt_id, chkpnt_iter, &
                     chkpnt_nac, chkpnt_new, chkpnt_rootener, chkpnt_rootener2, chkpnt_rootidx, chkpnt_rootmap, Iter_all
logical(kind=iwp) :: have_CI, have_NAC, have_RootEner2, have_RootMap
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
  call get_environment_variable('EMIL_InLoop',level)
  if ((level == '0') .or. (level == '1')) level = ''
  filename = basename//trim(level)

  if (((Iter > 1) .or. (IRC == -1)) .and. mh5_is_hdf5(filename)) then
    create = .false.
    chkpnt_id = mh5_open_file_rw(filename)
    chkpnt_iter = mh5_open_attr(chkpnt_id,'ITERATIONS')
    chkpnt_ener = mh5_open_dset(chkpnt_id,'ENERGIES')
    chkpnt_rootener = mh5_open_dset(chkpnt_id,'ROOT_ENERGIES')
    chkpnt_coor = mh5_open_dset(chkpnt_id,'COORDINATES')
    chkpnt_new = mh5_open_dset(chkpnt_id,'CENTER_COORDINATES')
    chkpnt_force = mh5_open_dset(chkpnt_id,'FORCES')
    chkpnt_hess = mh5_open_dset(chkpnt_id,'HESSIAN')
    ! conical intersection / MECI data: created only when the run has it (see Chkpnt_init), so a reopen must check for
    ! existence first -- mh5_open_dset does not, and would leave a poisoned handle that abends on first use.
    have_CI = mh5_exists_dset(chkpnt_id,'GRADIENT_DIFFERENCE')
    if (have_CI) then
      chkpnt_gd = mh5_open_dset(chkpnt_id,'GRADIENT_DIFFERENCE')
      chkpnt_rootidx = mh5_open_dset(chkpnt_id,'ROOT_INDICES')
    end if
    have_NAC = have_CI .and. mh5_exists_dset(chkpnt_id,'NAC')
    if (have_NAC) then
      chkpnt_nac = mh5_open_dset(chkpnt_id,'NAC')
      chkpnt_appnadc = mh5_open_dset(chkpnt_id,'APPROX_NADC')
    end if
    ! ROOT_MAPPING: created only for a Track run (see Chkpnt_init), so a reopen must check for existence first, exactly as
    ! for the CI/MECI datasets above.
    have_RootMap = mh5_exists_dset(chkpnt_id,'ROOT_MAPPING')
    if (have_RootMap) chkpnt_rootmap = mh5_open_dset(chkpnt_id,'ROOT_MAPPING')
    ! ROOT_ENERGIES_2: created only for a two-RunFile run (see Chkpnt_init), so a reopen must check for existence first,
    ! exactly as for the CI/MECI and ROOT_MAPPING datasets above.
    have_RootEner2 = mh5_exists_dset(chkpnt_id,'ROOT_ENERGIES_2')
    if (have_RootEner2) chkpnt_rootener2 = mh5_open_dset(chkpnt_id,'ROOT_ENERGIES_2')
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
  use Slapaf_Info, only: AtomLbl, Coor, dMass, dMEPStep, EDiffZero, iCoSet, iState, MEP, NADC, nDimBC, nStab, rMEP, Smmtrc, Track, &
                         TwoRunFiles
  use stdalloc, only: mma_allocate, mma_deallocate
  character :: lIrrep(24)
  integer(kind=iwp) :: Columbus, dsetid, i, j, k, mAtom, nRoots, nRoots2
  logical(kind=iwp) :: Found
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

  ! number of roots and per-root energies
  ! LDV: NROOTS is assumed fixed for the life of the file; unlike NSYM/NATOMS_UNIQUE in Chkpnt_open, it is not re-checked on
  ! reopen, since there is no mechanism by which the number of computed roots changes mid-optimization.
  nRoots = 1
  call Qpg_iScalar('Number of roots',Found)
  if (Found) call Get_iScalar('Number of roots',nRoots)
  call mh5_init_attr(chkpnt_id,'NROOTS',nRoots)
  chkpnt_rootener = mh5_create_dset_real(chkpnt_id,'ROOT_ENERGIES',2,[nRoots,0],dyn=.true.)
  call mh5_init_attr(chkpnt_rootener,'DESCRIPTION','Energies of all computed roots for all iterations, matrix of size '// &
                     '[ITERATIONS,NROOTS]; column order is assigned per geometry and a column need not hold the same state '// &
                     'across iterations -- ROOT_MAPPING resolves this when present, but its absence means the reordering was '// &
                     "not tracked, not that the order is fixed; on a TWO_RUNFILES run these are the active RunFile's roots "// &
                     "only -- RUNFILE2's roots are in ROOT_ENERGIES_2 when that dataset is present -- so ROOT_INDICES column "// &
                     '1 must not be used to index this dataset')

  ! root mapping: only meaningful (and only created) for a Track run, where the solver can reassign which root occupies
  ! which ROOT_ENERGIES column between iterations
  have_RootMap = Track
  if (have_RootMap) then
    chkpnt_rootmap = mh5_create_dset_int(chkpnt_id,'ROOT_MAPPING',2,[nRoots,0],dyn=.true.)
    call mh5_init_attr(chkpnt_rootmap,'DESCRIPTION',"Mapping from each root's original index at iteration 1 to its current "// &
                       'slot, matrix of size [ITERATIONS,NROOTS]: entry [k,i] is the slot occupied at iteration k by the root '// &
                       'that was index i at iteration 1')
  end if

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

  ! conical intersection / MECI data, created only when the run actually has it, following DESYM_FACTORS above rather than
  ! HESSIAN: an absent dataset says "not applicable", whereas one full of fill values cannot be told apart from a real result.
  ! Two conditions, deliberately different. have_CI is the wider one: Gx0 is a difference of two computed gradients on a
  ! same-spin CI run and on a two-RunFile crossing run alike, and on the latter it is the only two-state data in the file, so
  ! gating GRADIENT_DIFFERENCE on NADC would drop it exactly where it matters most. ROOT_INDICES and APPROX_NADC are
  ! per-iteration datasets rather than scalar attributes: ApproxNADC is reset every invocation and set only for that
  ! invocation's coupling failure, so an attribute would record the last geometry and label every earlier fallback as a real
  ! coupling; iState is RootMap-translated and can change across a root crossing.
  ! LDV: have_CI/have_NAC are decided once here (or on reopen, from dataset existence) and not re-verified every iteration;
  ! there is no mechanism by which iState(2) or NADC would flip mid-run for a file that already has these datasets, mirroring
  ! the NROOTS assumption above.
  have_CI = (iState(2) /= 0)
  if (have_CI) then
    call mh5_init_attr(chkpnt_id,'NADC',merge(1,0,NADC))
    call mh5_init_attr(chkpnt_id,'EDIFF_ZERO',merge(1,0,EDiffZero))
    call mh5_init_attr(chkpnt_id,'TWO_RUNFILES',merge(1,0,TwoRunFiles))

    ! per-root energies from RUNFILE2, present only for a two-RunFile run, where ROOT_ENERGIES alone cannot express the
    ! second wavefunction's roots
    have_RootEner2 = TwoRunFiles
    if (have_RootEner2) then
      nRoots2 = 1
      call NameRun('RUNFILE2')
      call Qpg_iScalar('Number of roots',Found)
      if (Found) call Get_iScalar('Number of roots',nRoots2)
      call NameRun('#Pop')
      call mh5_init_attr(chkpnt_id,'NROOTS_2',nRoots2)
      chkpnt_rootener2 = mh5_create_dset_real(chkpnt_id,'ROOT_ENERGIES_2',2,[nRoots2,0],dyn=.true.)
      call mh5_init_attr(chkpnt_rootener2,'DESCRIPTION','Energies of all computed roots of the wavefunction copied to RUNFILE2 '// &
                         'for a two-RunFile crossing run (e.g. the CASPT2 states of a CASPT2 gradient job, not necessarily '// &
                         'RASSCF), matrix of size [ITERATIONS,NROOTS_2]; ROOT_INDICES column 1 (Fortran index 2) selects '// &
                         'within this array, while ROOT_INDICES column 0 selects within ROOT_ENERGIES')
    end if

    chkpnt_gd = mh5_create_dset_real(chkpnt_id,'GRADIENT_DIFFERENCE',3,[3,size(Coor,2),0],dyn=.true.)
    call mh5_init_attr(chkpnt_gd,'DESCRIPTION','Gradient difference between the two states referenced by ROOT_INDICES (column '// &
                       '1 gradient minus column 0 gradient), matrix of size [ITERATIONS,NATOMS_UNIQUE,3], stored with '// &
                       'iteration varying slowest, then atom index')

    ! LDV: for a two-RunFile job iState is (active-RunFile root, RUNFILE2 root), not sorted -- see process_gradients.F90:
    ! L36 zeroes both, L64-75 sets iState(1) from the active RunFile, L104-119 sets iState(2) from RUNFILE2 with no min/max
    ! sort (unlike the same-spin CI case, where iState(1) ends up the higher root and iState(2) the lower).
    chkpnt_rootidx = mh5_create_dset_int(chkpnt_id,'ROOT_INDICES',2,[2,0],dyn=.true.)
    call mh5_init_attr(chkpnt_rootidx,'DESCRIPTION','State-pair indices for the two-state calculation, matrix of size '// &
                       '[ITERATIONS,2], holding Fortran 1-based root numbers -- a 0-based reader must subtract one; column 0 '// &
                       'is the higher root and column 1 the lower root, unless TWO_RUNFILES is set, in which case column 0 is '// &
                       "the active RunFile's root and column 1 is RUNFILE2's root, and the pair is not sorted; the NADC "// &
                       'attribute is 1 when a coupling derivative vector was computed for this pair, in which case NAC is '// &
                       'present, and 0 when it was not, i.e. a minimum-energy crossing point rather than a conical '// &
                       'intersection; EDIFF_ZERO is 1 when the energy difference is constrained to zero and 0 when a fixed '// &
                       'nonzero gap is sought')

    call Get_iScalar('Columbus',Columbus)
    ! Columbus /= 1 because in Columbus mode the block that fills NAC is skipped and the array keeps its zero fill; the flag
    ! is queried live where the array is filled, not taken from the input.
    have_NAC = NADC .and. (Columbus /= 1)
    if (have_NAC) then
      chkpnt_nac = mh5_create_dset_real(chkpnt_id,'NAC',3,[3,size(Coor,2),0],dyn=.true.)
      call mh5_init_attr(chkpnt_nac,'DESCRIPTION','Nonadiabatic coupling derivative vector between the two states referenced '// &
                         'by ROOT_INDICES, or (when APPROX_NADC is set for that iteration) a normalized dimensionless '// &
                         'branching-plane vector instead, matrix of size [ITERATIONS,NATOMS_UNIQUE,3], stored with iteration '// &
                         'varying slowest, then atom index; the NADC and EDIFF_ZERO attributes are described with ROOT_INDICES')

      chkpnt_appnadc = mh5_create_dset_int(chkpnt_id,'APPROX_NADC',1,[0],dyn=.true.)
      call mh5_init_attr(chkpnt_appnadc,'DESCRIPTION','Flag (0/1) per iteration: whether NAC for that iteration is an '// &
                         'approximate branching-plane vector rather than a true coupling derivative, vector of size [ITERATIONS]')
    end if
  else
    have_NAC = .false.
    have_RootEner2 = .false.
  end if

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
  use Slapaf_Info, only: ApproxNADC, Cx, Energy, Gx, Gx0, iState, iter, NAC, nDimBC, RootMap, TwoRunFiles
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp) :: i, ij, j, nRoots, nRoots2
  logical(kind=iwp) :: Found, FoundRoots
  real(kind=wp), allocatable :: Hss_X(:), RootEner(:), RootEner2(:)

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
  ! per-root energies
  nRoots = 1
  call Qpg_iScalar('Number of roots',FoundRoots)
  if (FoundRoots) call Get_iScalar('Number of roots',nRoots)
  call mma_allocate(RootEner,nRoots)
  call Get_dArray('Last energies',RootEner,nRoots)
  call mh5_resize_dset(chkpnt_rootener,[nRoots,Iter_all])
  call mh5_put_dset(chkpnt_rootener,RootEner,[nRoots,1],[0,Iter_all-1])
  call mma_deallocate(RootEner)
  ! per-root energies from RUNFILE2, on a two-RunFile run
  ! LDV: fCopy('RUNBACK','RUNFILE') at rlxctl.F90:331 restores the active RunFile from its saved snapshot but leaves
  ! RUNFILE2 untouched, so on an lNmHss run this row may come from a different geometry than ROOT_ENERGIES.
  if (TwoRunFiles) then
    call NameRun('RUNFILE2')
    nRoots2 = 1
    call Qpg_iScalar('Number of roots',Found)
    if (Found) call Get_iScalar('Number of roots',nRoots2)
    call mma_allocate(RootEner2,nRoots2)
    call Get_dArray('Last energies',RootEner2,nRoots2)
    call NameRun('#Pop')
    call mh5_resize_dset(chkpnt_rootener2,[nRoots2,Iter_all])
    call mh5_put_dset(chkpnt_rootener2,RootEner2,[nRoots2,1],[0,Iter_all-1])
    call mma_deallocate(RootEner2)
  end if
  ! root mapping
  if (have_RootMap) then
    call mh5_resize_dset(chkpnt_rootmap,[nRoots,Iter_all])
    call mh5_put_dset(chkpnt_rootmap,RootMap,[nRoots,1],[0,Iter_all-1])
  end if
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
  ! conical intersection / MECI data
  if (have_CI) then
    call mh5_resize_dset(chkpnt_gd,[3,size(Cx,2),Iter_all])
    call mh5_put_dset(chkpnt_gd,Gx0(:,:,Iter),[3,size(Cx,2),1],[0,0,Iter_all-1])
    call mh5_resize_dset(chkpnt_rootidx,[2,Iter_all])
    call mh5_put_dset(chkpnt_rootidx,iState,[2,1],[0,Iter_all-1])
    if (have_NAC) then
      call mh5_resize_dset(chkpnt_nac,[3,size(Cx,2),Iter_all])
      call mh5_put_dset(chkpnt_nac,NAC(:,:,Iter),[3,size(Cx,2),1],[0,0,Iter_all-1])
      call mh5_resize_dset(chkpnt_appnadc,[Iter_all])
      call mh5_put_dset(chkpnt_appnadc,[merge(1,0,ApproxNADC)],[1],[Iter_all-1])
    end if
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
