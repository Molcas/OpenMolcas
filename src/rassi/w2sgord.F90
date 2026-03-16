!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine W2SGORD(SGS,CIS,MWS2W,NLIST,KWALK,ICNUM)
! Purpose: Given a list of bit-packed total walks,
! translate this into a list of elements of a CI array.
! MXCPI= Max nr of case numbers packed in one integer.
! Dereference SGS:

use definitions, only: iwp
use gugx, only: SGStruct, CIStruct
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp) MWS2W(*), NLIST
integer(kind=iwp) KWALK(*), ICNUM(NLIST)
integer(kind=iwp), parameter :: MXCPI = 15
integer(kind=iwp), allocatable :: ICS(:)
integer(kind=iwp) nLev, nVert, MidLev, MVSta, nMidV, NIPWLK, MIPWLK

nLev = SGS%nLev
nVert = SGS%nVert
MidLev = SGS%MidLev
MVSta = SGS%MVSta
! Dereference CIS:

nMidV = CIS%nMidV
NIPWLK = CIS%nIpWlk
! Nr of integers used to store each total walk:
MIPWLK = 1+(NLEV-1)/MXCPI
! Allocate scratch space for case numbers:
call mma_allocate(ICS,NLEV,Label='ICS')
call W2SGORD1(NLEV,NVERT,NMIDV,NIPWLK,SGS%ISM,MIDLEV,MVSTA,CIS%IOCSF,CIS%NOW,CIS%IOW,SGS%DOWN,SGS%MAW,ICS,MWS2W,MIPWLK,NLIST, &
              KWALK,ICNUM)
call mma_deallocate(ICS)

end subroutine W2SGORD
