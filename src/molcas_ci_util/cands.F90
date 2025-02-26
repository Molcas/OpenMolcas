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

module CandS
! Definition of c and sigma
!
! iCSM  : Symmetry of CI state?
! iSSM  : Symmetry of Sigma state?
! iCSpc : Number of CI spaces
! iSSpc : Number of Sigma spaces

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: ICSM, ICSPC, ISSM, ISSPC

public :: ICSM, ICSPC, ISSM, ISSPC

end module CandS
