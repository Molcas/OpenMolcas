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

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
#include "periodic_table.fh"
character*2 Element(*)
real*8, allocatable :: Chrg(:)

call Get_iScalar('Unique atoms',nAtoms)
call mma_Allocate(Chrg,nAtoms)

call Get_dArray('Nuclear charge',Chrg,nAtoms)

do i=1,nAtoms
  iElement_Nr = int(Chrg(i))
  if ((iElement_Nr >= 0) .and. (iElement_Nr <= Num_Elem)) then
    Element(i) = PTab(iElement_Nr)
  else
    Element(i) = ' X'
  end if
end do
call mma_deallocate(Chrg)

return

end subroutine Get_Name
