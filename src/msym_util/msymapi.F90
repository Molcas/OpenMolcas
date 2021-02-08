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
! Copyright (C) 2015, Marcus Johansson                                 *
!***********************************************************************

subroutine fmsym_create_context(ctx)
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(out) :: ctx
integer(kind=iwp) :: ret
! INT cmsym_create_context(msym_context *pctx, int *err)
call cmsym_create_context(ctx,ret)
if (ret /= 0) then
  call WarningMessage(2,'Failed to create symmetry context')
  call Abend()
end if
return
end subroutine fmsym_create_context

!***********************************************************************

subroutine fmsym_set_elements(ctx)
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp) :: nAtoms
real(kind=wp), allocatable :: Coord(:)
call Get_nAtoms_All(nAtoms)
call mma_allocate(Coord,3*nAtoms,label='Coord')
call Get_Coord_All(Coord,nAtoms)
call fmsym_set_ele_orb(ctx,nAtoms,Coord)
call mma_deallocate(Coord)
return
end subroutine fmsym_set_elements

!***********************************************************************

subroutine fmsym_release_context(ctx)
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp) :: ret
! INT cmsym_release_context(msym_context *pctx, int*err)
call cmsym_release_context(ctx,ret)
return
end subroutine fmsym_release_context

!***********************************************************************

subroutine fmsym_set_ele_orb(ctx,nAtoms,Coord)
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: Coord(3,nAtoms)
#include "LenIn.fh"
character(len=LenIn), allocatable :: AtomLabel(:)
integer(kind=iwp), allocatable :: basis_ids(:), nBas(:)
integer(kind=iwp) :: nSym, nMO, nCMO, iSym, ret

call mma_allocate(AtomLabel,nAtoms,label='AtomLabel')
call Get_LblCnt_All(AtomLabel)
call Get_iScalar('nSym',nSym)
call mma_allocate(nBas,nSym,label='nBas')
call Get_iArray('nBas',nBas,nSym)

if (nSym /= 1) then
  call WarningMessage(2,'MSYM can only be used with group c1')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

nMO = 0
nCMO = 0
do iSym=1,nSym
  nMO = nMO + nBas(iSym)
  nCMO = nCMO + nBas(iSym)**2
end do
call mma_deallocate(nBas)

call mma_allocate(basis_ids,4*nMO,label='basis_ids')
call Get_iArray('Basis IDs',basis_ids,(4*nMO))
! INT cmsym_set_elements(msym_context *pctx, INT *pel, INT *puel, char *uelement, double xyz[][3], INT *paol, INT basis_ids[][4], int *err)
call cmsym_set_elements(ctx,nAtoms,(LENIN),AtomLabel,Coord,nMO,basis_ids,ret)
call mma_deallocate(basis_ids)
call mma_deallocate(AtomLabel)
call xflush(u6)
if (ret /= 0) then
  call WarningMessage(2,'Failed to set elements')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

return
end

!***********************************************************************

subroutine fmsym_find_symmetry(ctx)
use Definitions, only: wp, iwp, u6
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp) :: ret
character(len=6) :: PGName
! INT cmsym_find_symmetry(msym_context *pctx, char pgname[6], int *err)
call cmsym_find_symmetry(ctx,PGName,ret)
if (ret /= 0) then
  call WarningMessage(2,'Failed to find symmetry')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if
write(u6,*) 'Found Point Group:'
write(u6,*) PGName
return
end subroutine fmsym_find_symmetry

!***********************************************************************

Subroutine fmsym_symmetrize_molecule(ctx)
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp) :: lFN, ret
character(len=256) :: FN
call PrgmTranslate('MSYMOUT',FN,lFN)
FN = FN(1:lFN)//char(0)
! INT cmsym_symmetrize_molecule(msym_context *pctx, char *outfile, INT *err)
call cmsym_symmetrize_molecule(ctx,FN,ret)
if (ret /= 0) then
  call WarningMessage(2,'Failed to symmetrize molecule')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if
return
end subroutine fmsym_symmetrize_molecule

!***********************************************************************

subroutine fmsym_generate_orbital_subspaces(ctx)
#ifdef _HDF5_
use mh5, only: mh5_create_file, mh5_init_attr, mh5_create_dset_real, mh5_create_dset_int, mh5_create_dset_str, mh5_put_dset, &
               mh5_close_dset, mh5_close_file
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6
implicit none
real(kind=wp), intent(inout) :: ctx
integer(kind=iwp) :: iSym, nSym, nMO, nCMO, nIrr, LuOrb, Dummy(1), ret
integer(kind=iwp), allocatable :: nBas(:), IrrIds(:), IrrInd(:)
character(len=80) :: Title
character(len=8), allocatable :: irrep_strings(:)
real(kind=wp), allocatable :: CAO(:), Occ(:)
#ifdef _HDF5_
integer(kind=iwp) :: fileid, dsetid
#endif
integer(kind=iwp), external :: isfreeunit

call Get_iScalar('nSym',nSym)
call mma_allocate(nBas,nSym,label='nBas')
call Get_iArray('nBas',nBas,nSym)

nMO = 0
nCMO = 0
do iSym=1,nSym
  nMO = nMO + nBas(iSym)
  nCMO = nCMO + nBas(iSym)**2
end do

call mma_allocate(CAO,nCMO,label='CAO')
call mma_allocate(Occ,nMO,label='Occ')
call mma_allocate(IrrIds,nMO,label='IrrIds')
call mma_allocate(IrrInd,nMO,label='IrrInd')
call mma_allocate(irrep_strings,nMO,label='irrep_strings')

! INT cmsym_generate_orbital_subspaces(msym_context *pctx, INT *l, double c[*l][*l], INT irrep_ids[*l], INT irrep_ind[*l], INT *err)
call cmsym_generate_orbital_subspaces(ctx,nMO,CAO,IrrIds,IrrInd,nIrr,irrep_strings,ret)
write(u6,*) 'Irrep indices='
write(u6,'(5i3)') IrrInd(:)
write(u6,*) 'Irrep ids='
write(u6,'(5i3)') IrrIds(:)

if (ret /= 0) then
  call WarningMessage(2,'Failed to generate SALCs')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

Occ(:) = Zero

Title = 'Orbital Subspaces'
LuOrb = isfreeunit(50)
call WrVec('MSYMAORB',LuOrb,'CO',nSym,nBas,nBas,CAO,Occ,Dummy,IrrInd,Title)
#ifdef _HDF5_
fileid = mh5_create_file('MSYMH5')
call run2h5_molinfo(fileid)
call one2h5_ovlmat(fileid,nSym,nBas)
call one2h5_fckint(fileid,nSym,nBas)
! mocoef
dsetid = mh5_create_dset_real(fileid,'MO_VECTORS',1,[nCMO])
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Coefficients of the SALCs as produced by MSYM, arranged as blocks of size [NBAS(i)**2], i=1,#irreps')
call mh5_put_dset(dsetid,CAO)
call mh5_close_dset(dsetid)
! mooc
dsetid = mh5_create_dset_real(fileid,'MO_OCCUPATIONS',1,[nMO])
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Dummy occupation numbers arranged as blocks of size [NBAS(i)], i=1,#irreps')
call mh5_put_dset(dsetid,Occ)
call mh5_close_dset(dsetid)
! moene
dsetid = mh5_create_dset_real(fileid,'MO_ENERGIES',1,[nMO])
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Dummy orbital energies arranged as blocks of size [NBAS(i)], i=1,#irreps')
call mh5_put_dset(dsetid,Occ)
call mh5_close_dset(dsetid)
! supsym
dsetid = mh5_create_dset_int(fileid,'SUPSYM_IRREP_IDS',1,[nMO])
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Super-symmetry ids as produced by MSYM, arranged as blocks of size [NBAS(i)], i=1,#irreps')
call mh5_put_dset(dsetid,IrrIds)
call mh5_close_dset(dsetid)
dsetid = mh5_create_dset_int(fileid,'SUPSYM_IRREP_INDICES',1,[nMO])
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Super-symmetry indices as produced by MSYM, arranged as blocks of size [NBAS(i)], i=1,#irreps')
call mh5_put_dset(dsetid,IrrInd)
call mh5_close_dset(dsetid)
! irrep_labels
dsetid = mh5_create_dset_str(fileid,'SUPSYM_IRREP_LABELS',1,[nIrr],8)
call mh5_init_attr(dsetid, 'DESCRIPTION', &
                   'Super-symmetry labels as produced by MSYM, arranged as array of size i=1,#supsym_irreps')
call mh5_put_dset(dsetid,irrep_strings)
call mh5_close_dset(dsetid)
call mh5_close_file(fileid)
#endif

call mma_deallocate(nBas)
call mma_deallocate(CAO)
call mma_deallocate(Occ)
call mma_deallocate(IrrInd)
call mma_deallocate(IrrIds)
call mma_deallocate(irrep_strings)

return
end subroutine fmsym_generate_orbital_subspaces

!***********************************************************************

Subroutine fmsym_symmetrize_orbitals(ctx,CIO)
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(inout) :: ctx, CIO(*)
integer(kind=iwp) :: iSym, nSym, nMO, nCMO, ret
integer(kind=iwp), allocatable :: nBas(:)

call Get_iScalar('nSym',nSym)
call mma_allocate(nBas,nSym,label='nBas')
call Get_iArray('nBas',nBas,nSym)

if (nSym /= 1) then
  call WarningMessage(2,'MSYM can only be used with group c1')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

nMO = 0
nCMO = 0
do iSym=1,nSym
  nMO = nMO + nBas(iSym)
  nCMO = nCMO + nBas(iSym)**2
end do
call mma_deallocate(nBas)

! INT cmsym_symmetrize_orbitals(msym_context *pctx, INT *l, double c[*l][*l], INT *err)
call cmsym_symmetrize_orbitals(ctx,nMO,CIO,ret)
if (ret /= 0) then
  call WarningMessage(2,'Failed to symmetrize orbitals')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

return
end subroutine fmsym_symmetrize_orbitals

!***********************************************************************

Subroutine fmsym_symmetrize_orb_file(ctx,INPORB)
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
implicit none
real(kind=wp), intent(inout) :: ctx
character(len=*), intent(in) :: INPORB
integer(kind=iwp) :: iSym, iBas, nSym, nMO, nCMO, LuOrb, iWarn, iErr, i, j, ret
integer(kind=iwp), allocatable :: nBas(:), indType(:,:), TIND(:)
real(kind=wp), allocatable :: CIO(:), Occ(:), E(:)
character(len=80) :: Title
integer(kind=iwp), external :: isfreeunit

call Get_iScalar('nSym',nSym)
call mma_allocate(nBas,nSym,label='nBas')
call mma_allocate(indType,7,nSym,label='indType')
call Get_iArray('nBas',nBas,nSym)

if (nSym /= 1) then
  call WarningMessage(2,'MSYM can only be used with group c1')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if

nMO = 0
nCMO = 0
do iSym=1,nSym
   nMO = nMO + nBas(iSym)
   nCMO = nCMO + nBas(iSym)**2
end do

call mma_allocate(CIO,nCMO,label='nCMO')
call mma_allocate(Occ,nMO,label='Occ')
call mma_allocate(E,nMO,label='E')
call mma_allocate(TIND,nMO,label='TIND')

LuOrb = isfreeunit(50)
call RdVec(INPORB,LuOrb,'COEI',nSym,nBas,nBas,CIO,Occ,E,TIND,Title,iWarn,iErr)
! INT cmsym_symmetrize_orbitals(msym_context *pctx, INT *l, double c[*l][*l], INT *err)
call cmsym_symmetrize_orbitals(ctx,nMO,CIO,ret)

if (ret /= 0) then
  call WarningMessage(2,'Failed to symmetrize orbitals')
  ! INT cmsym_release_context(msym_context *pctx, int*err)
  call cmsym_release_context(ctx,ret)
  call Abend()
end if
indType(:,:) = 0
i = 1
do iSym=1,nSym
  do iBas=1,nBas(iSym)
    j = TIND(i)
    indType(j,iSym) = indType(j,iSym)+1
    i = i+1
  end do
end do
Title = 'Symmetrized Orbitals'
LuOrb = isfreeunit(50)
Call WrVec('MSYMMORB',LuOrb,'COEI',nSym,nBas,nBas,CIO,Occ,E,indType,Title)

call mma_deallocate(CIO)
call mma_deallocate(Occ)
call mma_deallocate(E)
call mma_deallocate(TIND)
call mma_deallocate(nBas)
call mma_deallocate(indType)

return
end subroutine fmsym_symmetrize_orb_file
