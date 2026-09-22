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
! PURPOSE: A SPIN-COUPLED CSF IS SPECIFIED BY NCL CLOSED SHELL
!          AND NOP OPEN SHELL AND OCCUPATION VECTORS ICL AND IOP,
!          RESPECTIVELY. THE SPIN COUPLING IS STORED IN THE
!          VECTOR ISPIN. TRANSLATE THESE DATA INTO THE
!          CORRESPONDING STEP VECTOR.
subroutine STEPVEC(ICL,IOP,NCL,NOP,ISPIN,NORB,IWALK)

use Definitions, only: iwp

implicit none

integer(iwp), intent(in)  :: NCL, NOP, NORB
integer(iwp), intent(in)  :: ICL(NCL), IOP(NOP)
integer(iwp), intent(in)  :: ISPIN(NOP)
integer(iwp), intent(out) :: IWALK(NORB)

integer(iwp) :: iClosed, iOpen
!
!  Start from an empty occupation pattern
!
IWALK(:) = 0
!
!  Mark doubly occupied orbitals
!
do iClosed=1,NCL
   IWALK(ICL(iClosed)) = 3
enddo
!
!  Mark singly occupied orbitals
!
do iOpen=1,NOP
   if (ISPIN(iOpen) == 1) then
      IWALK(IOP(iOpen)) = 1
   else
      IWALK(IOP(iOpen)) = 2
   endif
enddo

end subroutine STEPVEC
