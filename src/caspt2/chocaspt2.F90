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

module ChoCASPT2

use definitions, only: iwp

! Comment: numcho_pt2, InfVec_N2_PT2, and MaxVec_PT2 are copies
! of corresponding values in module Cholesky.
! Values are transferred at the beginning of caspt2.

implicit none

integer(kind=iwp) Lsplit(8), nIsplit(8), nAsplit(8), nksh(8), nkes(8), npsh(8), npes(8), numcho_pt2(8), iALGO, InfVec_N2_PT2, &
                  MaxVec_PT2, IF_CHO, NCHSPC, NHTSPC, NFTSPC, NFTSPC_TOT, MXNVC, MXCHARR

type ChoType
  integer(kind=iwp), allocatable :: Unit(:)
  integer(kind=iwp), allocatable :: ip(:)
  integer(kind=iwp), allocatable :: np(:)
  integer(kind=iwp), allocatable :: sp(:)
end type ChoType

type(ChoType) Stuff(8)

end module ChoCASPT2
