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

!#define _DEBUGPRINT_
subroutine PNT2DM(I12SM,NSMOB,IPSM,JPSM,IJSM,ISM2,IPNTR)
! Pointer to two dimensional array
!
! =====
! Input
! =====
! I12SM : /= 0 => restrict to lower half
!         == 0 => complete matrix
! NSMOB : Number of orbital symmetries
! IPSM  : Number of orbitals per symmetry for index 1
! JPSM  : Number of orbitals per symmetry for index 2
! IJSM  : Symmetry of two index array
!
! =======
! Output
! =======
! IPNTR : Pointer to block with first index of given symmetry
!         = 0 indicates forbidden block
! ISM2  : symmetry of second index for given first index

use Symmetry_Info, only: Mul
use Index_Functions, only: nTri_Elem
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: I12SM, NSMOB, IPSM(NSMOB), JPSM(NSMOB), IJSM
integer(kind=iwp), intent(out) :: ISM2(NSMOB), IPNTR(NSMOB)
integer(kind=iwp) :: IOFF, ISM, JSM

IPNTR(1:NSMOB) = 0
ISM2(1:NSMOB) = 0
IOFF = 1
do ISM=1,NSMOB
  JSM = Mul(ISM,IJSM)
  if (JSM == 0) exit
  if ((I12SM == 0) .or. (ISM >= JSM)) then
    ! Allowed block
    IPNTR(ISM) = IOFF
    ISM2(ISM) = JSM
    if ((I12SM > 0) .and. (ISM == JSM)) then
      IOFF = IOFF+nTri_Elem(IPSM(ISM))
    else
      IOFF = IOFF+IPSM(ISM)*JPSM(JSM)
    end if
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' dimension of two-dimensional array ',IOFF-1
write(u6,*) ' Pointer'
call IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
write(u6,*) ' Symmetry of other array'
call IWRTMA(ISM2,1,NSMOB,1,NSMOB)
#endif

end subroutine PNT2DM
