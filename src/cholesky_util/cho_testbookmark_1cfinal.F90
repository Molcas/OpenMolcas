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

subroutine Cho_TestBookmark_1CFinal()

use ChoArr, only: iAtomShl
use stdalloc, only: mma_deallocate

implicit none

if (allocated(iAtomShl)) call mma_deallocate(iAtomShl)

end subroutine Cho_TestBookmark_1CFinal
