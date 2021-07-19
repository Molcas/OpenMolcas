!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

function PageNo(iRoot)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the page number of a vector                              *
!                                                                      *
!     calling arguments:                                               *
!     iRoot   : integer                                                *
!               root number                                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use davctl_mod, only: istart, n_Roots, nvec
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: PageNo
integer(kind=iwp), intent(in) :: iRoot
integer(kind=iwp) :: itmp1
#include "rasdim.fh"

itmp1 = iRoot
if (iRoot > n_Roots) then
  itmp1 = n_Roots+mod(istart+iRoot-n_Roots-1,nvec-n_Roots)+1
end if
PageNo = itmp1

return

end function PageNo
