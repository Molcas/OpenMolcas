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

subroutine MSTOW1(NSYM,NLEV,NVERT,NMIDV,NIPWLK,NWALK,MIDLEV,ICS,NOW,IOW,IWALK,IUP,IDOWN,MAW,MWS2W)
! Purpose: From the list of packed up- and downwalks, construct
! the table MWS2W, such that MAW sums can be translated to the
! corresponding walks of the Split-GUGA.

use definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NSYM, NLEV, NVERT, NMIDV, NIPWLK, NWALK, MIDLEV
integer(kind=iwp), intent(out) :: ICS(NLEV)
integer(kind=iwp), intent(in) :: NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
integer(kind=iwp), intent(in) :: IWALK(NIPWLK*NWALK)
integer(kind=iwp), intent(in) :: IDOWN(NVERT,0:3), IUP(NVERT,0:3)
integer(kind=iwp), intent(in) :: MAW(NVERT,0:3)
integer(kind=iwp), intent(out) :: MWS2W(NWALK)
integer(kind=iwp) MV, ISYUP, NUP, IUOFF, IUW, IUWTOT, MS, IUV, LEV, IC, ISYDWN, NDWN, IDOFF, IDW, IDWTOT, IDV

do MV=1,NMIDV
  do ISYUP=1,NSYM
    NUP = NOW(1,ISYUP,MV)
    if (NUP == 0) cycle
    IUOFF = IOW(1,ISYUP,MV)/NIPWLK
    do IUW=1,NUP
      IUWTOT = IUOFF+IUW
      ! Unpack upper walk to ICS()
      call UPKWLK(NLEV-MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IUWTOT-1)),ICS(MIDLEV+1))
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
      call UPKWLK(MIDLEV,NIPWLK,1,IWALK(1+NIPWLK*(IDWTOT-1)),ICS)
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
