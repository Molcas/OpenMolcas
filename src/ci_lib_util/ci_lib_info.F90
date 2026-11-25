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
!> (CI_Lib_States), for the calling programs: a copy made at the
!> initialization. The programs read it here and never in the CI
!> environments (SGUGA, Faroald), which keep their own data.
module CI_Lib_Info

use Molcas, only: MxLev, MxSym
use CI_Space_Def, only: CI_Space, CI_Space_Input
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: CI_Lib_MxSpace = 4
integer(kind=iwp) :: iq

!> nSym, nActEl, iSpin: symmetry, electrons, multiplicity of the CI space;
!> nLev: number of levels (active orbitals) of the split graph, MidLev its
!> midlevel; IFRAS: RAS flag of the graph (> 2: more than one RAS
!> partition is occupied); nCSF(iSym): CSFs per total symmetry;
!> Level(iOrb): level of active orbital iOrb, L2Act(iLev): active orbital
!> of level iLev, ISm(iLev): symmetry of level iLev.
type CI_Lib_Space_Info
  logical(kind=iwp) :: Active = .false.
  integer(kind=iwp) :: nSym = 0, nActEl = 0, iSpin = 0, nLev = 0, MidLev = 0, IFRAS = 0, nRsPrt = 0
  integer(kind=iwp) :: nCSF(MxSym) = 0
  ! (Level and L2Act default to the identity, as in the split graph: the
  ! programs may read them before the setup, e.g. CASPT2's pt2wfn_init)
  integer(kind=iwp) :: Level(MxLev) = [(iq,iq=1,MxLev)], L2Act(MxLev) = [(iq,iq=1,MxLev)], ISm(MxLev) = 0
  ! GAS: the state holding the relaxed space of this one (0: none), see
  ! CI_SG_Rel_Space
  integer(kind=iwp) :: iRel = 0
  ! nDet(iSym): dimension of the determinant space of the CI back end per
  ! symmetry (Lucia: determinants or spin combinations; without Lucia
  ! nCSF); nCSF_Kept: CSFs of the CI symmetry kept by HEXS/DEXS
  integer(kind=iwp) :: nDet(MxSym) = 0, nCSF_Kept = 0
  ! nSD(iSym): Slater determinants per symmetry (Lucia only, for the
  ! output); Spin_Comb: nDet counts spin combinations (Lucia's IDC = 2)
  integer(kind=iwp) :: nSD(MxSym) = 0
  logical(kind=iwp) :: Spin_Comb = .false.
  ! STSym: the symmetry of the CI state (CI_Setup_SG; one graph covers all
  ! symmetries)
  integer(kind=iwp) :: STSym = 0
end type CI_Lib_Space_Info

type(CI_Lib_Space_Info) :: CI_Lib_Space(CI_Lib_MxSpace)

!> GAS: the state index of the relaxed space of the RASSCF CI space: the
!> GAS bounds widened by one electron at every partition boundary, so that
!> it holds every single excitation of a GAS CSF and products of
!> excitation operators (H, 2-RDM) are exact there for GAS vectors.
integer(kind=iwp), parameter :: CI_SG_Rel_Space = CI_Lib_MxSpace

public :: CI_Lib_Space, CI_Lib_Space_Info, CI_Lib_MxSpace, CI_SG_Rel_Space

!> The CI space of RASSCF in the common description of all back ends
!> (CI_Space_Def: partitions with cumulative bounds, removed HEXS/DEXS
!> configurations), built by the CI interface (CI_Space_Setup).
type(CI_Space) :: Space_CI
!> The HEXS/DEXS requests of the input (set by RASSCF's input processing),
!> from which CI_Space_Setup builds the removals of Space_CI.
type(CI_Space_Input) :: Space_Inp
public :: Space_CI, Space_Inp

end module CI_Lib_Info
