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

subroutine NYPART(iExtra,nPart,COORD,rStart,nCent,iSeed)

use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "warnings.h"
integer(kind=iwp) :: iExtra, nPart, nCent, iSeed
real(kind=wp) :: COORD(MxCen*MxPut,3), rStart
integer(kind=iwp) :: I, IA, IIN, iMAXVAR, IND, IVARV
real(kind=wp) :: DR, DX, DY, DZ, R2, RLIM, RLIM2, rlm, X, Y, Z
real(kind=wp), external :: Ranf

! Preparing some numbers
iMAXVAR = 100*iExtra**2
IIN = 1
IVARV = 0
RLIM = rStart-8.5_wp
rlm = rlim*Two
RLIM2 = RLIM**Two
! Check if not all molecules been put out in reasonable time
outer: do
  if (IVARV >= iMAXVAR) then
    write(u6,*) 'Failure to add particles. Try to increase the dielectric radie or change the random seed.'
    call Quit(_RC_GENERAL_ERROR_)
  end if
  ! Here is the random change relative the first user defined water
  DX = ranf(iseed)*RLM-rlim
  DY = ranf(iseed)*RLM-rlim
  DZ = ranf(iseed)*RLM-rlim
  DR = DX*DX+DY*DY+DZ*DZ
  IVARV = IVARV+1
  if (DR > RLIM2) cycle outer
  IND = IIN+NPART
  X = DX+COORD(1,1)
  Y = DY+COORD(1,2)
  Z = DZ+COORD(1,3)
  ! Check so that two water molecules do not come too close...
  do I=1,IND*NCENT,NCENT
    R2 = (COORD(I,1)-X)**2+(COORD(I,2)-Y)**2+(COORD(I,3)-Z)**2
    if (R2 < 60.0_wp) cycle outer
  end do
  ! ...and if they do not then shove its coordinates in variable
  IIN = IIN+1
  IA = (IIN+NPART-1)*NCENT
  do I=1,NCENT
    COORD(IA+I,1) = COORD(I,1)+DX
    COORD(IA+I,2) = COORD(I,2)+DY
    COORD(IA+I,3) = COORD(I,3)+DZ
  end do
  ! Check if all water molecules have been put where they should
  if (IIN >= IEXTRA) exit outer
end do outer
! Give nPart its new value and exit gracefully
NPART = NPART+IEXTRA

return

end subroutine NYPART
