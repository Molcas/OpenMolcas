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

subroutine SG_MSTOW(SGS,CIS,MWS2W,nSym)

use sguga, only: CIStruct, SGStruct
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

contains

subroutine MSTOW1(NSYM,NLEV,NVERT,NMIDV,NIPWLK,NWALK,MIDLEV,ICS,NOW,IOW,IWALK,IUP,IDOWN,MAW,MWS2W)
! Purpose: From the list of packed up- and downwalks, construct
! the table MWS2W, such that MAW sums can be translated to the
! corresponding walks of the Split-GUGA.

implicit none
integer(kind=iwp), intent(in) :: NSYM, NLEV, NVERT, NMIDV, NIPWLK, NWALK, MIDLEV, NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV), &
                                 IWALK(NIPWLK*NWALK), IUP(NVERT,0:3), IDOWN(NVERT,0:3), MAW(NVERT,0:3)
integer(kind=iwp), intent(out) :: ICS(NLEV), MWS2W(NWALK)
integer(kind=iwp) :: IC, IDOFF, IDV, IDW, IDWTOT, ISYDWN, ISYUP, IUOFF, IUV, IUW, IUWTOT, LEV, MS, MV, NDWN, NUP

do MV=1,NMIDV
  do ISYUP=1,NSYM
    NUP = NOW(1,ISYUP,MV)
    if (NUP == 0) cycle
    IUOFF = IOW(1,ISYUP,MV)/NIPWLK
    do IUW=1,NUP
      IUWTOT = IUOFF+IUW
      ! Unpack upper walk to ICS()
      call SG_UPKWLK(NLEV-MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IUWTOT-1)),ICS(MIDLEV+1))
      MS = 0
      IUV = 1
      do LEV=NLEV,MIDLEV+1,-1
        IC = ICS(LEV)
        MS = MS+MAW(IUV,IC)
        IUV = IDOWN(IUV,IC)
      end do
      MWS2W(MS) = IUWTOT
    end do
  end do
end do

do MV=1,NMIDV
  do ISYDWN=1,NSYM
    NDWN = NOW(2,ISYDWN,MV)
    if (NDWN == 0) cycle
    IDOFF = IOW(2,ISYDWN,MV)/NIPWLK
    do IDW=1,NDWN
      IDWTOT = IDOFF+IDW
      ! Unpack lower walk to ICS()
      call SG_UPKWLK(MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IDWTOT-1)),ICS)
      MS = 0
      IDV = NVERT
      do LEV=1,MIDLEV
        IC = ICS(LEV)
        IUV = IUP(IDV,IC)
        MS = MS+MAW(IUV,IC)
        IDV = IUV
      end do
      MWS2W(MS) = IDWTOT
    end do
  end do
end do

end subroutine MSTOW1

end subroutine SG_MSTOW
