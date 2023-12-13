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
integer(kind=iwp), intent(in) :: iExtra, nCent
integer(kind=iwp), intent(inout) :: nPart, iSeed
real(kind=wp), intent(inout) :: COORD(3,(nPart+iExtra)*nCent)
real(kind=wp), intent(in) :: rStart
integer(kind=iwp) :: I, IA, IIN, iMAXVAR, IND, IVARV
real(kind=wp) :: DR, DX, DY, DZ, R2, RLIM, RLIM2, rlm, X, Y, Z
real(kind=wp), external :: Random_Molcas
#include "warnings.h"

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
  DX = Random_Molcas(iseed)*RLM-rlim
  DY = Random_Molcas(iseed)*RLM-rlim
  DZ = Random_Molcas(iseed)*RLM-rlim
  DR = DX*DX+DY*DY+DZ*DZ
  IVARV = IVARV+1
  if (DR > RLIM2) cycle outer
  IND = IIN+NPART
  X = DX+COORD(1,1)
  Y = DY+COORD(2,1)
  Z = DZ+COORD(3,1)
  ! Check so that two water molecules do not come too close...
  do I=1,IND*NCENT,NCENT
    R2 = (COORD(1,I)-X)**2+(COORD(2,I)-Y)**2+(COORD(3,I)-Z)**2
    if (R2 < 60.0_wp) cycle outer
  end do
  ! ...and if they do not then shove its coordinates in variable
  IIN = IIN+1
  IA = (IIN+NPART-1)*NCENT
  COORD(1,IA+1:NCENT) = COORD(1,1:NCENT)+DX
  COORD(2,IA+1:NCENT) = COORD(2,1:NCENT)+DY
  COORD(3,IA+1:NCENT) = COORD(3,1:NCENT)+DZ
  ! Check if all water molecules have been put where they should
  if (IIN >= IEXTRA) exit outer
end do outer
! Give nPart its new value and exit gracefully
NPART = NPART+IEXTRA

return

end subroutine NYPART
