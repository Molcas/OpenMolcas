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

module gas_data

use Molcas, only: MxGas, MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IGSOCCX(mxGAS,2) = 0, NGAS = 0, NGSSH(mxGAS,mxSym) = 0
logical(kind=iwp) :: iDoGas = .false.

public :: iDoGas, IGSOCCX, NGAS, NGSSH

end module gas_data
