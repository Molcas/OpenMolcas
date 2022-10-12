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

subroutine Get_Name(Element)

use Isotopes, only: PTab
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
character(len=2), intent(_OUT_) :: Element(*)
integer(kind=iwp) :: i, iElement_Nr, nAtoms
real(kind=wp), allocatable :: Chrg(:)

call Get_iScalar('Unique atoms',nAtoms)
call mma_Allocate(Chrg,nAtoms)

call Get_dArray('Nuclear charge',Chrg,nAtoms)

do i=1,nAtoms
  iElement_Nr = int(Chrg(i))
  if ((iElement_Nr >= lbound(PTab,1)) .and. (iElement_Nr <= ubound(PTab,1))) then
    Element(i) = PTab(iElement_Nr)
  else
    Element(i) = ' X'
  end if
end do
call mma_deallocate(Chrg)

return

end subroutine Get_Name
