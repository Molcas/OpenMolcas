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

module RASDim

use Definitions, only: iwp

implicit none
private

! Parameter definitions
! Note, these parameters define the size of the problems which
! can be treated by the RASSCF program. Changing any of these
! requires, that the program is totaly recompiled and linked.
!
! mxRef  max number of reference configurations in root selectioning
! mxIter max number of macro iterations
! mxCiIt max number of micro iterations for the CI section
! mxSxIt max number of micro iterations for the SX section
! mxTit  max number of title lines

integer(kind=iwp), parameter :: mxCiIt = 502, mxIter = 200, mxRef = 5, mxSxIt = 100, mxTit = 1

public :: mxCiIt, mxIter, mxRef, mxSxIt, mxTit

end module RASDim
