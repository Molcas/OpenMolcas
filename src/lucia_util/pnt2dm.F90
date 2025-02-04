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

subroutine PNT2DM(I12SM,NSMOB,NSMSX,OSXO,IPSM,JPSM,IJSM,ISM2,IPNTR,MXPOBS)
! Pointer to two dimensional array
!
! =====
! Input
! =====
! I12SM  : ne.0 => restrict to lower half
!          eq.0 => complete matrix
! NSMOB : Number of orbital symmetries
! NSMSX : Number of SX      symmetries
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

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
integer OSXO(MXPOBS,2*MXPOBS), IPSM(*), JPSM(*)
! Output
dimension IPNTR(*), ISM2(*)

call ISETVC(IPNTR,0,NSMOB)
call ISETVC(ISM2,0,NSMOB)
IOFF = 1
do ISM=1,NSMOB
  JSM = OSXO(ISM,IJSM)
  if (JSM == 0) goto 100
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
100 continue
end do

NTEST = 0
if (NTEST >= 1) write(u6,*) ' dimension of two-dimensional array ',IOFF-1
if (NTEST >= 5) then
  write(u6,*) ' Pointer'
  call IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
  write(u6,*) ' Symmetry of other array'
  call IWRTMA(ISM2,1,NSMOB,1,NSMOB)
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(NSMSX)

end subroutine PNT2DM
