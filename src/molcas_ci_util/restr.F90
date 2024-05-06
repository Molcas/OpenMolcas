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

subroutine RESTR(SGS)
! PURPOSE: PUT THE RAS CONSTRAINT TO THE DRT TABLE BY CREATING A MASK

use gugx, only: SGStruct
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, IV, IVD, IVV, LEV, MASK, N
integer(kind=iwp), parameter :: i_and(0:3,0:3) = reshape([0,0,0,0,0,1,0,1,0,0,2,2,0,1,2,3],[4,4]), &
                                i_or(0:3,0:3) = reshape([0,1,2,3,1,1,3,3,2,3,2,3,3,3,3,3],[4,4]), &
                                LTAB = 1, NTAB = 2

call mma_allocate(SGS%Ver,SGS%nVert0,Label='V11')

! LOOP OVER ALL VERTICES AND CHECK ON RAS CONDITIONS
! CREATE MASK

do IV=1,SGS%nVert0
  LEV = SGS%DRT0(IV,LTAB)
  N = SGS%DRT0(IV,NTAB)
  SGS%Ver(IV) = 0
  if ((LEV == SGS%LV1RAS) .and. (N >= SGS%LM1RAS)) SGS%Ver(IV) = 1
  if ((LEV == SGS%LV3RAS) .and. (N >= SGS%LM3RAS)) SGS%Ver(IV) = SGS%Ver(IV)+2
end do

! NOW LOOP FORWARDS, MARKING THOSE VERTICES CONNECTED FROM ABOVE.
! SINCE VER WAS INITIALIZED TO ZERO, NO CHECKING IS NEEDED.

do IV=1,SGS%nVert0-1
  IVV = SGS%Ver(IV)
  do IC=0,3
    ID = SGS%Down0(IV,IC)
    if (ID /= 0) SGS%Ver(ID) = i_or(SGS%Ver(ID),IVV)
  end do
end do

! THEN LOOP BACKWARDS. SAME RULES, EXCEPT THAT CONNECTIVITY
! SHOULD BE PROPAGATED ONLY ABOVE THE RESTRICTION LEVELS.

do IV=SGS%nVert0-1,1,-1
  LEV = SGS%DRT0(IV,LTAB)
  MASK = 0
  if (LEV > SGS%LV1RAS) MASK = 1
  if (LEV > SGS%LV3RAS) MASK = MASK+2
  IVV = SGS%Ver(IV)
  do IC=0,3
    ID = SGS%Down0(IV,IC)
    if (ID /= 0) then
      IVD = SGS%Ver(ID)
      IVV = i_or(IVV,i_and(MASK,IVD))
    end if
  end do
  SGS%Ver(IV) = IVV
end do

! WE ARE NOW INTERESTED ONLY IN VERTICES CONNECTED BOTH TO
! ALLOWED VERTICES FOR RAS-SPACE 1 AND RAS-SPACE 3.
! THOSE ARE NUMBERED IN ASCENDING ORDER, THE REST ARE ZEROED.

SGS%nVert = 0
do IV=1,SGS%nVert0
  if (SGS%Ver(IV) == 3) then
    SGS%nVert = SGS%nVert+1
    SGS%Ver(IV) = SGS%nVert
  else
    SGS%Ver(IV) = 0
  end if
end do
if (SGS%nVert == 0) call SysAbendMsg('Restr','No configuration was found\n','Check NACTEL, RAS1, RAS2, RAS3 values')
#ifdef _DEBUGPRINT_
write(u6,*) 'RESTR:'
write(u6,*) 'LV1RAS, LV3RAS, LM1RAS, LM3RAS=',SGS%LV1RAS,SGS%LV3RAS,SGS%LM1RAS,SGS%LM3RAS
do IV=1,SGS%nVert0
  write(u6,*) 'VER(:)=',SGS%Ver(IV)
end do
#endif

end subroutine RESTR
