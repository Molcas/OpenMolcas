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

subroutine GETSTEPVECTOR(NOW,IOW,MV,IDWN,IUP,ICS,nLev,nMidV)

use gugx, only: CIS, SGS
use general_data, only: NSYM
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nMidV, NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV), nLev
integer(kind=iwp), intent(inout) :: MV, IDWN, IUP
integer(kind=iwp), intent(out) :: ICS(NLEV)
integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICUP, ICUPOS, IDW0, IUW0, LEV, NDWN, NNN, NUP

! RECONSTRUCT THE CASE LIST

! ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
! WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.

NUP = NOW(1,1,MV)
NDWN = NOW(2,1,MV)
IUW0 = 1-CIS%nIpWlk+IOW(1,1,MV)
IDW0 = 1-CIS%nIpWlk+IOW(2,1,MV)
! determine the stepvector
ICDPOS = IDW0+IDWN*CIS%nIpWlk
ICDWN = CIS%ICASE(ICDPOS)
! unpack lower walk
NNN = 0
do LEV=1,SGS%MIDLEV
  NNN = NNN+1
  if (NNN == 16) then
    NNN = 1
    ICDPOS = ICDPOS+1
    ICDWN = CIS%ICASE(ICDPOS)
  end if
  IC1 = ICDWN/4
  ICS(LEV) = ICDWN-4*IC1
  ICDWN = IC1
end do
ICUPOS = IUW0+CIS%nIpWlk*IUP
ICUP = CIS%ICASE(ICUPOS)
! unpack upper walk
NNN = 0
do LEV=SGS%MIDLEV+1,NLEV
  NNN = NNN+1
  if (NNN == 16) then
    NNN = 1
    ICUPOS = ICUPOS+1
    ICUP = CIS%ICASE(ICUPOS)
  end if
  IC1 = ICUP/4
  ICS(LEV) = ICUP-4*IC1
  ICUP = IC1
end do

! compute the next set of indices
if (IUP == NUP) then
  if (IDWN == NDWN) then
    if (MV == NMIDV) then
      MV = 0
    else
      MV = MV+1
    end if
    IDWN = 1
  else
    IDWN = IDWN+1
  end if
  IUP = 1
else
  IUP = IUP+1
end if

end subroutine GETSTEPVECTOR
