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

subroutine PNT2DM(I12SM,NSMOB,OSXO,IPSM,JPSM,IJSM,ISM2,IPNTR,MXPOBS)
! Pointer to two dimensional array
!
! =====
! Input
! =====
! I12SM  : ne.0 => restrict to lower half
!          eq.0 => complete matrix
! NSMOB : Number of orbital symmetries
! OSXO  : Symmetry of orbital, SX => symmetry of other orbital
! IPSM : Number of orbitals per symmetry for index 1
! JPSM : Number of orbitals per symmetry for index 2
! IJSM  : Symmetry of two index array
!
! =======
! Output
! =======
! IPNTR : Pointer to block with first index of given symmetry
!         = 0 indicates forbidden block
! ISM2  : symmetry of second index for given first index

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: I12SM, NSMOB, MXPOBS, OSXO(MXPOBS,2*MXPOBS), IPSM(*), JPSM(*), IJSM, ISM2(*), IPNTR(*)
integer(kind=iwp) :: IOFF, ISM, JSM, NTEST

IPNTR(1:NSMOB) = 0
ISM2(1:NSMOB) = 0
IOFF = 1
do ISM=1,NSMOB
  JSM = OSXO(ISM,IJSM)
  if (JSM == 0) exit
  if ((I12SM == 0) .or. (ISM >= JSM)) then
    ! Allowed block
    IPNTR(ISM) = IOFF
    ISM2(ISM) = JSM
    if ((I12SM > 0) .and. (ISM == JSM)) then
      IOFF = IOFF+IPSM(ISM)*(IPSM(ISM)+1)/2
    else
      IOFF = IOFF+IPSM(ISM)*JPSM(JSM)
    end if
  end if
end do

NTEST = 0
if (NTEST >= 1) write(u6,*) ' dimension of two-dimensional array ',IOFF-1
if (NTEST >= 5) then
  write(u6,*) ' Pointer'
  call IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
  write(u6,*) ' Symmetry of other array'
  call IWRTMA(ISM2,1,NSMOB,1,NSMOB)
end if

end subroutine PNT2DM
