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

subroutine Get_Name_Full(Element)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
character(len=2), intent(_OUT_) :: Element(*)
#include "Molcas.fh"
integer(kind=iwp) :: i, nAtMM, nAtom
logical(kind=iwp) :: Found
character(len=LenIn), allocatable :: LabMM(:)

call Get_nAtoms_All(nAtom)
call Get_Name_All(Element)

call Qpg_cArray('MMO Labels',Found,nAtMM)
if (Found) then
  nAtMM = nAtMM/LENIN
  call mma_allocate(LabMM,nAtMM,label='MMO Labels')
  call Get_cArray('MMO Labels',LabMM,LENIN*nAtMM)
  do i=1,nAtMM
    Element(nAtom+i) = LabMM(i)(1:2)
    if (Element(nAtom+i)(2:2) == '_') Element(nAtom+i)(2:2) = ' '
  end do
  call mma_deallocate(LabMM)
end if

return

end subroutine Get_Name_Full
