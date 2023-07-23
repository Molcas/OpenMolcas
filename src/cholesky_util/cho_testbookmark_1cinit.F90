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

subroutine Cho_TestBookmark_1CInit(AllocatedHere)

use ChoArr, only: iAtomShl
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
logical(kind=iwp) :: AllocatedHere
#include "cholesky.fh"
integer(kind=iwp) :: irc

if (.not. allocated(iAtomShl)) then
  call mma_allocate(iAtomShl,nShell,Label='iAtomShl')
  irc = -1
  call Cho_SetAtomShl(irc,iAtomShl,size(iAtomShl))
  if (irc /= 0) then
    write(u6,'(A,I4)') 'Cho_TestBookmark_1Cinit: Cho_SetAtomShl returned',irc
    call Cho_Quit('shell-to-atom init failed!',104)
  end if
  AllocatedHere = .true.
else
  AllocatedHere = .false.
end if

end subroutine Cho_TestBookmark_1CInit
