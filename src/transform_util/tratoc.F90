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

module TraToc

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: nTraBuf = 9600, nTraToc = 106
integer(kind=iwp) :: iTraToc(nTraToc)

public :: iTraToc, nTraBuf, nTraToc

end module TraToc
