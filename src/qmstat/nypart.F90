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

subroutine NYPART(iExtra,nPart,COORD,Rstart,nCent,iSeed)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "warnings.h"
dimension COORD(MxCen*MxPut,3)
external Ranf

! Preparing some numbers
iMAXVAR = 100*iExtra**2
IN = 1
IVARV = 0
RLIM = RSTART-8.5d0
rlm = rlim*2.
RLIM2 = RLIM**2
! Check if not all molecules been put out in reasonable time
100 if (IVARV >= iMAXVAR) then
  write(6,*) 'Failure to add particles. Try to increase the dielectric radie or change the random seed.'
  call Quit(_RC_GENERAL_ERROR_)
end if
! Here is the random change relative the first user definied water
DX = ranf(iseed)*RLM-rlim
DY = ranf(iseed)*RLM-rlim
DZ = ranf(iseed)*RLM-rlim
DR = DX*DX+DY*DY+DZ*DZ
IVARV = IVARV+1
if (DR > RLIM2) GO TO 100
IND = IN+NPART
X = DX+COORD(1,1)
Y = DY+COORD(1,2)
Z = DZ+COORD(1,3)
! Check so that two water molecules do not come too close...
do I=1,IND*NCENT,NCENT
  R2 = (COORD(I,1)-X)**2+(COORD(I,2)-Y)**2+(COORD(I,3)-Z)**2
  if (R2 < 60d0) Go To 100
end do
! ...and if they do not then shove its coordinates in variable
IN = IN+1
IA = (IN+NPART-1)*NCENT
do I=1,NCENT
  COORD(IA+I,1) = COORD(I,1)+DX
  COORD(IA+I,2) = COORD(I,2)+DY
  COORD(IA+I,3) = COORD(I,3)+DZ
end do
! Check if all water molecules have been put where they should
if (IN < (IEXTRA)) GO TO 100
! Give nPart its new value and exit gracefully
NPART = NPART+IEXTRA

return

end subroutine NYPART
