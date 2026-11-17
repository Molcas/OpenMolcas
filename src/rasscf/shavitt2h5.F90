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

!> Write the split graph (Shavitt graph) of the CI space of state iState
!> to the group SHAVITT of the HDF5 file fileid, for Tools/Shavitt.
!>
!> The graph: the distinct-row table DRT (L, N, a, b, c per vertex; vertex 1
!> is the head at level nLev), the down-chain table DOWN (the vertex below
!> for the steps d = 0..3, 0: no arc) and the symmetry of each level.
!> The CSF order of CI_VECTORS: the counts and offsets of the half walks
!> (NOW, IOW) and of the CSF blocks (NOCSF, IOCSF), and the packed step
!> vectors of all half walks (ICASE, NPACK steps of 2 bits per integer,
!> NIPWLK integers per half walk, lowest level in the lowest bits).
!> Nothing is written if the graph has not been built (no SGUGA).
subroutine Shavitt2h5(fileid,iState)

#ifdef _HDF5_
use mh5, only: mh5_close_group, mh5_create_group, mh5_init_attr, mh5_init_dset
use SGUGA, only: CIS, nPack, SGS
use sxci, only: IDXSX
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: fileid, iState
#ifdef _HDF5_
integer(kind=iwp) :: grpid, nLev, nMidV, nRsPrt, nSym, nVert

if (.not. (allocated(SGS(iState)%DRT) .and. allocated(SGS(iState)%Down) .and. allocated(CIS(iState)%ICase))) return

nSym = SGS(iState)%nSym
nLev = SGS(iState)%nLev
nVert = SGS(iState)%nVert
nMidV = CIS(iState)%nMidV
nRsPrt = SGS(iState)%nRsPrt

grpid = mh5_create_group(fileid,'SHAVITT')

call mh5_init_attr(grpid,'DESCRIPTION','Split graph (Shavitt graph) of the CSF space and the CSF order of CI_VECTORS')
call mh5_init_attr(grpid,'NSYM',nSym)
call mh5_init_attr(grpid,'NACTEL',SGS(iState)%nActEl)
call mh5_init_attr(grpid,'ISPIN',SGS(iState)%iSpin)
call mh5_init_attr(grpid,'NLEV',nLev)
call mh5_init_attr(grpid,'NVERT',nVert)
call mh5_init_attr(grpid,'MIDLEV',SGS(iState)%MidLev)
call mh5_init_attr(grpid,'MVSTA',SGS(iState)%MVSta)
call mh5_init_attr(grpid,'MVEND',SGS(iState)%MVEnd)
call mh5_init_attr(grpid,'NMIDV',nMidV)
call mh5_init_attr(grpid,'IFRAS',SGS(iState)%IFRAS)
call mh5_init_attr(grpid,'NRSPRT',nRsPrt)
call mh5_init_attr(grpid,'NPACK',nPack)
call mh5_init_attr(grpid,'NIPWLK',CIS(iState)%nIpWlk)
call mh5_init_attr(grpid,'NWALK',CIS(iState)%nWalk)
call mh5_init_attr(grpid,'NUW',CIS(iState)%NUW)

call mh5_init_dset(grpid,'DRT',2,[nVert,5],pack(SGS(iState)%DRT,.true.))
call Descr('DRT','Distinct-row table, for each vertex: level L, electrons N and Paldus a, b, c; vertex 1 is the head')
call mh5_init_dset(grpid,'DOWN',2,[nVert,4],pack(SGS(iState)%Down,.true.))
call Descr('DOWN','Down-chain table, for each vertex the vertex below for the steps d=0,1,2,3 (0: no arc)')
call mh5_init_dset(grpid,'LEVEL_SYM',1,[nLev],SGS(iState)%ISm)
call Descr('LEVEL_SYM','Irrep of the orbital of each level 1..NLEV')
! RASSCF numbers the levels RAS/GAS partition by partition; IDXSX maps them to the active orbitals
call mh5_init_dset(grpid,'LEVEL_TO_ACT',1,[nLev],IDXSX(1:nLev))
call Descr('LEVEL_TO_ACT','Active orbital index of each level 1..NLEV')
call mh5_init_dset(grpid,'NRAS',2,[nSym,nRsPrt],pack(SGS(iState)%nRas(1:nSym,1:nRsPrt),.true.))
call Descr('NRAS','Number of orbitals of each irrep in each RAS/GAS partition')
call mh5_init_dset(grpid,'NRASEL',1,[nRsPrt],SGS(iState)%nRasEl)
call Descr('NRASEL','Minimum number of electrons in the partitions 1..k')
call mh5_init_dset(grpid,'NRASMX',1,[nRsPrt],SGS(iState)%nRasMx)
call Descr('NRASMX','Maximum number of electrons in the partitions 1..k (-1: no bound)')
call mh5_init_dset(grpid,'NCSF',1,[nSym],CIS(iState)%NCSF)
call Descr('NCSF','Number of CSFs of each total irrep')
call mh5_init_dset(grpid,'NOW',3,[2,nSym,nMidV],pack(CIS(iState)%NOW,.true.))
call Descr('NOW','Number of upper (1) and lower (2) half walks by irrep and midvertex')
call mh5_init_dset(grpid,'IOW',3,[2,nSym,nMidV],pack(CIS(iState)%IOW,.true.))
call Descr('IOW','Offsets in ICASE of the NOW blocks')
call mh5_init_dset(grpid,'NOCSF',3,[nSym,nMidV,nSym],pack(CIS(iState)%NOCSF,.true.))
call Descr('NOCSF','Number of CSFs by upper-walk irrep, midvertex and total irrep')
call mh5_init_dset(grpid,'IOCSF',3,[nSym,nMidV,nSym],pack(CIS(iState)%IOCSF,.true.))
call Descr('IOCSF','Offsets of the NOCSF blocks in the CI vector')
call mh5_init_dset(grpid,'ICASE',1,[CIS(iState)%nWalk*CIS(iState)%nIpWlk],CIS(iState)%ICase)
call Descr('ICASE','Packed step vectors of all half walks, NIPWLK integers per walk, NPACK steps of 2 bits per integer')

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
