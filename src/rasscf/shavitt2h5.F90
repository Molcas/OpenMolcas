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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> Write the definition of the split graph (Shavitt graph) of the CI space
!> of state iState to the group SHAVITT of the HDF5 file fileid, for
!> Tools/Shavitt.
!>
!> Together with NSYM, NACTEL, SPINMULT, LSYM and L2ACT of the file this
!> defines the graph and the CSF order of CI_VECTORS: Tools/Shavitt rebuilds
!> them as MKSGUGA does. Stored are the irrep of each level, the orbitals of
!> each irrep in the partitions (RAS/GAS spaces, levels numbered partition by
!> partition) with the cumulative electron bounds, and the midlevel; the
!> number of vertices and the CSFs per irrep let the tool check its graph.
!> Nothing is written without a graph (no SGUGA, or no active orbitals).
subroutine Shavitt2h5(fileid,iState)

#ifdef _HDF5_
use mh5, only: mh5_close_group, mh5_create_group, mh5_init_attr, mh5_init_dset
use SGUGA, only: CIS, SGS
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: fileid, iState
#ifdef _HDF5_
integer(kind=iwp) :: grpid, nLev, nRsPrt, nSym

if (.not. (allocated(SGS(iState)%DRT) .and. allocated(SGS(iState)%ISm) .and. allocated(CIS(iState)%NCSF))) return

nSym = SGS(iState)%nSym
nLev = SGS(iState)%nLev
nRsPrt = SGS(iState)%nRsPrt
if (nLev < 1) return  ! no active orbitals: no graph

grpid = mh5_create_group(fileid,'SHAVITT')

call mh5_init_attr(grpid,'DESCRIPTION','Definition of the split graph (Shavitt graph) of the CSF space; with NSYM, '// &
                   'NACTEL, SPINMULT, LSYM and L2ACT it gives the graph and the CSF order of CI_VECTORS (SGUGA, MKSGUGA)')
call mh5_init_attr(grpid,'MIDLEVEL',SGS(iState)%MidLev)
call mh5_init_attr(grpid,'NUM_VERTICES',SGS(iState)%nVert)

call mh5_init_dset(grpid,'LEVEL_IRREPS',1,[nLev],SGS(iState)%ISm(1:nLev))
call Descr('LEVEL_IRREPS','Irrep of the orbital of each level 1..nLev (level 1 at the bottom of the graph)')
call mh5_init_dset(grpid,'PARTITION_ORBITALS',2,[nSym,nRsPrt],pack(SGS(iState)%nRas(1:nSym,1:nRsPrt),.true.))
call Descr('PARTITION_ORBITALS','Number of orbitals of each irrep in each partition (RAS/GAS space)')
call mh5_init_dset(grpid,'PARTITION_MIN_ELECTRONS',1,[nRsPrt],SGS(iState)%nRasEl(1:nRsPrt))
call Descr('PARTITION_MIN_ELECTRONS','Minimum number of electrons in the partitions 1..k together')
call mh5_init_dset(grpid,'PARTITION_MAX_ELECTRONS',1,[nRsPrt],SGS(iState)%nRasMx(1:nRsPrt))
call Descr('PARTITION_MAX_ELECTRONS','Maximum number of electrons in the partitions 1..k together (-1: no bound)')
call mh5_init_dset(grpid,'CSFS_PER_IRREP',1,[nSym],CIS(iState)%NCSF(1:nSym))
call Descr('CSFS_PER_IRREP','Number of CSFs of each total irrep')

call mh5_close_group(grpid)

contains

subroutine Descr(dset,text)
  use mh5, only: mh5_close_dset, mh5_open_dset, mh5_init_attr
  character(len=*), intent(in) :: dset, text
  integer(kind=iwp) :: dsetid
  dsetid = mh5_open_dset(grpid,dset)
  call mh5_init_attr(dsetid,'DESCRIPTION',text)
  call mh5_close_dset(dsetid)
end subroutine Descr

#else
#include "macros.fh"
unused_var(fileid)
unused_var(iState)
#endif

end subroutine Shavitt2h5
