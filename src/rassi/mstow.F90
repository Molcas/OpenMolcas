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

subroutine MSTOW(SGS,CIS,MWS2W,nSym)

use gugx, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(_OUT_) :: MWS2W(*)
integer(kind=iwp), intent(in) :: nSym
integer(kind=iwp) :: MIDLEV, NIPWLK, NLEV, NMIDV, NVERT, NWALK
integer(kind=iwp), allocatable :: ICS(:)

NLEV = SGS%nLev
NVERT = SGS%nVert
MIDLEV = SGS%MidLev

NMIDV = CIS%nMidV
NIPWLK = CIS%nIpWlk
NWALK = CIS%nWalk
call mma_allocate(ICS,NLEV,Label='ICS')
call MSTOW1(NSYM,NLEV,NVERT,NMIDV,NIPWLK,NWALK,MIDLEV,ICS,CIS%NOW,CIS%IOW,CIS%ICase,SGS%UP,SGS%DOWN,SGS%MAW,MWS2W)
call mma_deallocate(ICS)

end subroutine MSTOW
