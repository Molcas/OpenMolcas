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

subroutine SetPrLev(IPRGLB_IN,IPRLOC_IN)

use PrintLevel, only: DEBUG, SILENT, USUAL
use output_ras, only: IPRGLB, IPRLOC
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPRGLB_IN, IPRLOC_IN(7)
integer(kind=iwp) :: I
logical(kind=iwp), external :: REDUCE_PRT

! The local print levels are the maximum of the requested global and
! local ones, except that if any of IPRGLB or IPRLOC(I) is zero
! (meaning silence), then IPRLOC(I) is set to zero.
IPRGLB = IPRGLB_IN
IPRLOC(:) = 0
if (IPRGLB_IN /= 0) then
  do I=1,7
    if (IPRLOC_IN(I) > 0) IPRLOC(I) = max(IPRGLB_IN,IPRLOC_IN(I))
  end do
end if
! If inside an optimization loop, set down the print
! level unless we *really* want a lot of output.
if (REDUCE_PRT()) then
  IPRGLB = max(IPRGLB-USUAL,SILENT)
  IPRLOC(:) = max(IPRLOC(:)-USUAL,SILENT)
end if

if (IPRLOC(1) >= DEBUG) then
  write(u6,*) ' SetPrLev: Print levels have been set to'
  write(u6,*) '  Global print level IPRGLB=',IPRGLB
  write(u6,*) '  Individual sections print levels, IPRLOC:'
  write(u6,'(1x,7I5)') (IPRLOC(I),I=1,7)
end if

end subroutine SetPrLev
