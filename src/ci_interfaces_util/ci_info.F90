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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The CI information of the states set up through the CI interface
!> (CI_States), for the calling programs: a copy made at the
!> initialization. The programs read it here and never in the CI
!> environments (SGUGA, Faroald), which keep their own data.
module CI_Info

use Molcas, only: MxLev, MxSym
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MxCIState = 3

!> nSym, nActEl, iSpin: symmetry, electrons, multiplicity of the CI space;
!> nLev: number of levels (active orbitals) of the split graph, MidLev its
!> midlevel; IFRAS: RAS flag of the graph (> 2: more than one RAS
!> partition is occupied); nCSF(iSym): CSFs per total symmetry;
!> Level(iOrb): level of active orbital iOrb, L2Act(iLev): active orbital
!> of level iLev, ISm(iLev): symmetry of level iLev.
type CI_State_Info
  logical(kind=iwp) :: Active = .false.
  integer(kind=iwp) :: nSym = 0, nActEl = 0, iSpin = 0, nLev = 0, MidLev = 0, IFRAS = 0, nRsPrt = 0
  integer(kind=iwp) :: nCSF(MxSym) = 0
  integer(kind=iwp) :: Level(MxLev) = 0, L2Act(MxLev) = 0, ISm(MxLev) = 0
end type CI_State_Info

type(CI_State_Info) :: CI_St(MxCIState)

public :: CI_St, CI_State_Info, MxCIState

end module CI_Info
