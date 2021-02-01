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

module Surfacehop_globals

use Definitions, only: iwp

implicit none
private

logical(kind=iwp) :: lH5Restart
character(len=180) :: File_H5Res

public :: lH5Restart, File_H5Res

end module Surfacehop_globals
