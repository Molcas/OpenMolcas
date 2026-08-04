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
! Purpose: From the list of packed up- and downwalks, construct
! the table MWS2W, such that MAW sums can be translated to the
! corresponding walks of the Split-GUGA.

use sguga, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(_OUT_) :: MWS2W(*)
integer(kind=iwp), intent(in) :: nSym
integer(kind=iwp) :: IC, IDOFF, IDV, IDW, IDWTOT, ISYDWN, ISYUP, IUOFF, IUV, IUW, IUWTOT, LEV, MS, MV, NDWN, NUP
integer(kind=iwp), allocatable :: ICS(:)

call mma_allocate(ICS,SGS%nLev,Label='ICS')

do MV=1,CIS%nMidV
  do ISYUP=1,NSYM
    NUP = CIS%NOW(1,ISYUP,MV)
    if (NUP == 0) cycle
    IUOFF = CIS%IOW(1,ISYUP,MV)/CIS%nIpWlk
    do IUW=1,NUP
      IUWTOT = IUOFF+IUW
      ! Unpack upper walk to ICS()
      call SG_UPKWLK(SGS%nLev-SGS%MidLev,CIS%nIpWlk,1,CIS%ICase(1+CIS%nIpWlk*(IUWTOT-1)),ICS(SGS%MidLev+1))
      MS = 0
      IUV = 1
      do LEV=SGS%nLev,SGS%MidLev+1,-1
        IC = ICS(LEV)
        MS = MS+SGS%MAW(IUV,IC)
        IUV = SGS%DOWN(IUV,IC)
      end do
      MWS2W(MS) = IUWTOT
    end do
  end do
end do

do MV=1,CIS%nMidV
  do ISYDWN=1,NSYM
    NDWN = CIS%NOW(2,ISYDWN,MV)
    if (NDWN == 0) cycle
    IDOFF = CIS%IOW(2,ISYDWN,MV)/CIS%nIpWlk
    do IDW=1,NDWN
      IDWTOT = IDOFF+IDW
      ! Unpack lower walk to ICS()
      call SG_UPKWLK(SGS%MidLev,CIS%nIpWlk,1,CIS%ICase(1+CIS%nIpWlk*(IDWTOT-1)),ICS)
      MS = 0
      IDV = SGS%nVert
      do LEV=1,SGS%MidLev
        IC = ICS(LEV)
        IUV = SGS%UP(IDV,IC)
        MS = MS+SGS%MAW(IUV,IC)
        IDV = IUV
      end do
      MWS2W(MS) = IDWTOT
    end do
  end do
end do

call mma_deallocate(ICS)

end subroutine SG_MSTOW
