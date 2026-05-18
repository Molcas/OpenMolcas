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

! Comment: numcho_pt2, InfVec_N2_PT2, and MaxVec_PT2 are copies
! of corresponding values in module Cholesky.
! Values are transferred at the beginning of caspt2.

use Definitions, only: iwp

implicit none
private

type ChoType
  integer(kind=iwp), allocatable :: Unt(:)
  integer(kind=iwp), allocatable :: ip(:)
  integer(kind=iwp), allocatable :: np(:)
  integer(kind=iwp), allocatable :: sp(:)
end type ChoType

integer(kind=iwp) :: iALGO, MaxVec_PT2, MXCHARR, MXNVC, nAsplit(8), NCHSPC, NFTSPC, NFTSPC_TOT, NHTSPC, nIsplit(8), nksh(8), &
                     npsh(8), numcho_pt2(8)
type(ChoType) :: Stuff(8)

public :: iALGO, MaxVec_PT2, MXCHARR, MXNVC, NASplit, NCHSPC, NFTSPC, NFTSPC_TOT, NHTSPC, NISplit, nksh, npsh, NumCho_PT2, Stuff

end module ChoCASPT2
